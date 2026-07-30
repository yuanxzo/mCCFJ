classdef mCCFJ
% A MATLAB Package for calculating seismic ambient noise cross-correlation and frequency-Bessel transformation.
%   
% In this package, these are four main functions: mCCFJ.correlate, mCCFJ.transform, mCCFJ.inversion and mCCFJ.attenuation.
% mCCFJ.correlate is used to calculate the cross-correlation or
%       cross-coherency function of seismic waveforms. To ensure the
%       efficiency of calculation, the function is calculated in the
%       frequency domain, and the cross-correlation between any two
%       stations in the same time window is calculated in the way of matrix
%       parallelism. If necessary, GPU acceleration can be used.
% mCCFJ.transform is used for dispersion analysis of cross-correlation
%       function. We provide a variety of frequency-wavenumber domain
%       transformation methods to deal with different data, which is up to
%       you. Compared with previous methods, we provide an enhanced version
%       of frequency-Bessel transform here, i.e. spatial windowed
%       frequency-Bessel transform. This new method can make the energy of
%       the dispersion spectrum more concentrated and reduce spatial
%       artifacts, which is beneficial to the analysis of seismic wave
%       phase velocity and attenuation. This feature can be used by
%       specifying ops.win='hamming_1', 'hamming_2' or 'hamming_half' in
%       mCCFJ.transform. In order to ensure the efficiency of computing,
%       GPU acceleration can be used when necessary.
% mCCFJ.inversion is the inverse transformation program of mCCFJ.transform,
%       which allows for the inverse transformation of some of the
%       dispersion energy in the frequency-velocity domain back to obtain
%       the cross-correlation signal in the frequency-space domain. This
%       ability helps with denoising, mode separation, etc.
% mCCFJ.attenuation extracts surface‑wave attenuation by analyzing the
%       amplitudes of the cross-coherency outputs from mCCFJ.inversion
%       using an improved coherence‑fitting method.
% 
% For more details, refer to the following paper, and thank you for quoting
% it if the mCCFJ program brings convenience to your research. Reference:
%       Yang, B., Meng, H., Yuan, S., & Chen, X. (2025). 
%       Reliable Multimodal Attenuation Estimation of Surface Waves Using Diffuse Ambient Noise: Theory and Applications.
%       Journal of Geophysical Research: Solid Earth, 130, e2025JB031418. https://doi.org/10.1029/2025JB031418
%
% To_begin_with = 'mCCFJ.Help'
% Timestamp: 2026-07-15
% Copyright © 2026 Bo Yang (杨博). All rights reserved.

methods (Static)

    % Correlate seismic waveforms
    function CC=correlate(wave,rloc,Fs,ops)
        arguments
            wave (:,:) {mustBeNumeric}
            rloc (:,:) {mustBeNumeric} = [1:1:length(wave(1,:));zeros(1,length(wave(1,:)));zeros(1,length(wave(1,:)))]
            Fs (1,1) {mustBeNumeric,mustBePositive} = 1
            ops.AX {mustBeMember(ops.AX,["latlon","xyz"])} = "latlon"
            ops.CL (1,1) {mustBeInteger,mustBePositive} = length(wave(:,1))
            ops.OL (1,1) {mustBeInteger,mustBeNonnegative} = 0
            ops.FL (1,1) {mustBeInteger,mustBeNonnegative} = 0
            ops.NT {mustBeMember(ops.NT,["No","OneBit"])} = "No"
            ops.NF {mustBeMember(ops.NF,["No","PSD","ABS"])} = "No"
            ops.FM (1,1) {mustBeNumeric,mustBePositive} = Fs/2
            ops.FD {mustBeInteger,mustBePositive}=1    % Frequency point dilution
            ops.TP {mustBeMember(ops.TP,["No","Hann","tukeywin_5","tukeywin_10"])} = "tukeywin_5"
            ops.RR (:,6) = []
            ops.GPU{mustBeMember(ops.GPU,["No","Yes"])} = "No"
        end
        CC=struct;
    
        % 解码可选参数
        [npts,nsta]=size(wave,[1 2]);
        if nsta<=1
            error("At least two seismic waveform data are required!")
        end
        if npts<=2
            error("The number of sampling points in the waveform data is less than 2!")
        end
        if ops.CL > npts
            error("The length of the sliding window for correlation calculation exceeds the total length of the waveform, ending calculation!");
        end
        if ops.OL >= ops.CL
            error("The overlap length of the sliding window exceeds the length of the sliding window, ending calculation!");
        end
        if isempty(rloc)
            rloc=zeros(2,nsta);
            rloc(1,:)=1:nsta;
            rloc(2,:)=0;
            ops.AX = "xyz";
        else
            [m,n]=size(rloc,[1,2]);
            if m<2 || n~=nsta
                error("The size of the station location data does not match the number of seismic waveform data!");
            end
        end
        if ops.FM > Fs/2
            warning("The maximum frequency of interest is greater than half of the sampling rate, and now 'FM' is set to 'Fs/2'");
            ops.FM = Fs/2;
        end
        if ops.FL==0            
            ops.FL=ops.CL;
            % ops.FL=2^nextpow2(ops.CL);
        end
        if ops.FL<ops.CL
            error("'FL' must be larger than or equal to 'CL'!");
        end
        if ops.FD~=1
            warning("The spectra of waveforms will be decimated at equal intervals to speed up computation. However, this operation causes time‑domain aliasing, making the resulting time‑domain cross‑correlation unreliable.")
        end
    
        % 确定互相关滑动窗的起始点的位置及个数
        index = 1 : ops.CL-ops.OL : npts-ops.CL+1;
        if npts - index(end) + 1 < ops.CL
            index(end)=[];
        end
        num_of_win = length(index);
        len_of_win = ops.FL;
    
        % 计算所有可能的台站对的站间距离
        if isempty(ops.RR)
            RR=mCCFJ.distances(rloc,ops.AX); % RR(:,1) 为距离, RR(:,2) 为第一个台站的索引, RR(:,3) 为第二个台站的索引, RR(:,4) 为方位角
        else
            if isa(ops.RR,"cell")
                RR=ops.RR;
            else
                error("options.RR type mismatch, please calculate it according to 'mCCFJ.distances.'")
            end
        end
        num_of_ccr = length(RR(:,1));
        st1=int64(cell2mat(RR(:,2))); 
        st2=int64(cell2mat(RR(:,3)));

        % One-bit
        if strcmp(ops.NT,"OneBit")
            wave = sign(wave);
        end
    
        % 频率序列
        Freq = 0:Fs/(ops.FL):Fs/2;
        finx = Freq<=ops.FM;
        freq = Freq(finx);
        freq = freq(1:ops.FD:end);
        
        % 分配内存
        acf = zeros(length(freq),nsta);
        ccf = randn(length(freq),num_of_ccr);ccf = ccf + 1i.*ccf;
        astack = zeros(1, nsta);        % 自相关叠加次数
        cstack = zeros(1, num_of_ccr);  % 互相关叠加次数
        
        % taper
        switch ops.TP
            case "No"
                taper = 1;
            case "Hann"
                taper = hann(ops.CL);
            case "tukeywin_5"
                taper = tukeywin(ops.CL,0.05);
            case "tukeywin_10"
                taper = tukeywin(ops.CL,0.10);
        end

        % 检查GPU是否可用
        if canUseGPU() && strcmp(ops.GPU,"Yes")
            cb=whos("ccf",'wave');
            gm=gpuDevice().AvailableMemory;
            if 2*cb(1).bytes+cb(2).bytes < gm*0.9 
                wave=gpuArray(wave);
                acf=gpuArray(acf);
                ccf=gpuArray(ccf);
                st1=gpuArray(st1);
                st2=gpuArray(st2);
                index=gpuArray(index); 
                len_of_win=gpuArray(len_of_win);
                num_of_win=gpuArray(num_of_win);
                taper=gpuArray(taper);
                finx=gpuArray(finx);
                astack=gpuArray(astack);
                cstack=gpuArray(cstack);

                GPUtype=1;
            elseif 2*cb(1).bytes < gm*0.8 
                acf=gpuArray(acf);
                ccf=gpuArray(ccf);
                st1=gpuArray(st1);
                st2=gpuArray(st2);
                index=gpuArray(index); 
                len_of_win=gpuArray(len_of_win);
                num_of_win=gpuArray(num_of_win);
                taper=gpuArray(taper);
                finx=gpuArray(finx);
                astack=gpuArray(astack);
                cstack=gpuArray(cstack); 

                GPUtype=1;
            else
                GPUtype=0;
                disp("Attempted to use the GPU, but graphics memory is insufficient. Run on CPU")
            end
        else
            GPUtype=0;
        end

        
        ccf(:,:)=0;
        % 逐个滑动窗计算互相关, 平均为最终结果
        for i=1:num_of_win
            wave_seg = wave(index(i):index(i)+ops.CL-1,:);

            % 如果有些道的数据含有0过多则不用此道
            temp = sum(sign(abs(wave_seg)),1);
            inx0 = temp<ops.CL*0.5;
            wave_seg(:,inx0)=0;

            % 如果仅至多一道波形有数据, 互相关无意义, 跳过计算
            temp = sign(sum(abs(wave_seg),1));
            if length(find(temp==1)) <= 1
                continue
            end

            if GPUtype==1
                wave_seg=gpuArray(wave_seg);
            end

            % 计算频谱
            wave_seg = fft(wave_seg.*taper,len_of_win,1)/len_of_win;
            wave_seg = wave_seg(1:floor(len_of_win/2)+1,:);
            wave_seg(2:end-1,:) = 2*wave_seg(2:end-1,:);

            % 只保留关心的最高频率以下的信息, 加快计算速度
            wave_seg = wave_seg(finx,:);

            % 如果只想关注频率域的结果，且为了加速计算，可考虑对频谱做等间距的抽稀处理
            wave_seg = wave_seg(1:ops.FD:end,:);

            % 计算自相关运算
            temp = real(conj(wave_seg).*wave_seg);
    
            % 叠加自相关
            inx1 = find(sum(abs(temp),1)>0);
            acf(:,inx1) = acf(:,inx1) + temp(:,inx1);
            astack(inx1) = astack(inx1) + 1;

            % 互相关 和 归一化
            cc_temp = conj(wave_seg(:,st1)).*wave_seg(:,st2);
            switch ops.NF
                case "No"   % 直接互相关
                    cc_norm = 1;
                case "PSD"  % PSD归一化
                    cc_norm = mean(temp(:,inx1),2); % 当前时窗的PSD
                    cc_norm = cc_norm + mean(abs(cc_norm),"all").*1e-6;
                case "ABS"  % 绝对值归一化，完全的谱白化
                    cc_norm = sqrt(temp(:,st1)).*sqrt(temp(:,st2));
                    cc_norm = cc_norm + mean(abs(cc_norm),"all").*1e-6;
            end
            cc_temp = cc_temp./cc_norm;
            
            % 标记各列数据的有效性
            inx2 = find(sum(abs(cc_temp),1)>0);
    
            % 叠加互相关
            ccf(:,inx2)  = ccf(:,inx2)  + cc_temp(:,inx2);
            cstack(inx2) = cstack(inx2) + 1;
        end
        
        % 所有滑动窗的结果取平均
        inx3 = find(astack>0);
        acf(:,inx3) = acf(:,inx3) ./ astack(inx3);
        inx3 = find(cstack>0);
        ccf(:,inx3) = ccf(:,inx3) ./ cstack(inx3);


        % 如果使用了FM，补齐到该有的长度
        acf(end+1:length(Freq),:)=0;
        ccf(end+1:length(Freq),:)=0;

        % 不符合要求的强行置为0
        inx4 = astack==0; acf(:,inx4) = 0;
        inx4 = cstack==0; ccf(:,inx4) = 0;
        
        % 计算时间域互相关
        cct = fftshift( myifft_cc(ccf,len_of_win), 1);

        % 输出结果    
        CC.acf = gather(acf);
        CC.ccf = gather(ccf);
        CC.cct = gather(cct);
        CC.ccr = RR;
        CC.freq = gather(Freq');
        CC.time=((1:length(CC.cct(:,1)))-(floor(length(CC.cct(:,1))/2) + 1))./Fs/ops.FD;        
        CC.info.astack = gather(astack);
        CC.info.cstack = gather(cstack);
        CC.info.num_of_win = gather(num_of_win);
        CC.info.len_of_win = gather(len_of_win);
    end

    % 计算互相关距离
    function RR=distances(rloc,AX,ops)
        arguments
            rloc (:,:) {mustBeNumeric}
            AX   {mustBeMember(AX,["latlon","xyz"])}
            ops.sta_name (1,:) = 1:length(rloc(1,:)) % 台站名默认为它们在rloc中的位置
        end
        if length(rloc(:,1))==2
            rloc(3,:)=0;
        end
        nsta=length(rloc(1,:));
        RR=cell((nsta*(nsta-1)/2),6);

        nccr = 0;
        if strcmp(AX,"xyz")     % 直角坐标
            for i = 1:nsta-1
                s1=rloc(:,i);
                s2=rloc(:,i+1:nsta);
                RR(nccr+1:nccr+nsta-i,1)=num2cell(sqrt(sum((s1-s2).^2,1)));
                RR(nccr+1:nccr+nsta-i,2)={i};
                RR(nccr+1:nccr+nsta-i,3)=num2cell(i+1:nsta);
                RR(nccr+1:nccr+nsta-i,4)=num2cell(mod(450-rad2deg(atan2(s2(2,:)-s1(2), s2(1,:)-s1(1))), 360));
                if isa(ops.sta_name,"double")
                    RR(nccr+1:nccr+nsta-i,5)={ops.sta_name(i)};
                    RR(nccr+1:nccr+nsta-i,6)=num2cell(ops.sta_name(i+1:nsta));
                elseif isa(ops.sta_name,"cell")
                    RR(nccr+1:nccr+nsta-i,5)=ops.sta_name(i);
                    RR(nccr+1:nccr+nsta-i,6)=ops.sta_name(i+1:nsta);
                elseif isa(ops.sta_name,"string")
                    RR(nccr+1:nccr+nsta-i,5)={ops.sta_name(i)};
                    for j=i+1:nsta
                        RR(nccr+j-i,6)={ops.sta_name(j)};
                    end
                end
                nccr=nccr+nsta-i;
            end            
        elseif strcmp(AX,"latlon") % 经纬度坐标
            for i = 1:nsta-1
                s1=rloc(:,i);
                for j = i+1:nsta
                    s2=rloc(:,j);
                    nccr=nccr+1;
                    [dist,  az]=distance(s1(2), s1(1), s2(2), s2(1), referenceEllipsoid('WGS 84'));
                    RR(nccr,1)={sqrt(dist.^2+(s2(3)-s1(3)).^2)};
                    RR(nccr,2)={i};
                    RR(nccr,3)={j};
                    RR(nccr,4)={az};      % 前一个台作观察后一个台的方位角
                    if isa(ops.sta_name,"double")
                        RR(nccr,5)={ops.sta_name(i)};
                        RR(nccr,6)={ops.sta_name(j)};
                    elseif isa(ops.sta_name,"cell")
                        RR(nccr,5)=ops.sta_name(i);
                        RR(nccr,6)=ops.sta_name(j);
                    elseif isa(ops.sta_name,"string")
                        RR(nccr,5)={ops.sta_name(i)};
                        RR(nccr,6)={ops.sta_name(j)};
                    end
                end
            end
        end

        [~, order] = sort(cell2mat(RR(:,1)));
        RR=RR(order,:);
    end


    % 频率-贝塞尔变换
    function FJ=transform(CC, c_range, f_bound, ops)
        arguments
            CC (1,1) struct
            c_range (1,:) {mustBeNumeric,mustBeNonzero}
            f_bound (1,2) {mustBeNumeric,mustBeNonnegative}
            ops.Fun {mustBeMember(ops.Fun,["J0","J1","H1","H2","H2_r","FK"])}="H2"
            ops.Win {mustBeMember(ops.Win,["No","Hamming_1","Hamming_2","Hamming_half","Gauss"])}="Hamming_1"
            ops.GPU {mustBeMember(ops.GPU,["No","Yes"])}="No"
            ops.Num (1,1){mustBeNumeric,mustBeNonnegative}=0;      
        end
        FJ=struct;
        
        if strcmp(ops.GPU,"Yes") && ops.Num>0
            warning('CPU + GPU parallel execution may slow down the program. It is recommended to enable only one mode.')
        end

        % CC 是 correlate 函数计算得出的, 如果不是，按下面形式合成CC也可
        frq = CC.freq;       % 频率域互相关依赖的频率序列
        if isa(CC.ccr,"cell")
            ccr = cell2mat(CC.ccr(:,1)); % 互相关对的距离
        elseif isa(CC.ccr,"double")
            ccr = CC.ccr(:,1);
        end
        ccr = ccr(:);
        ccf = CC.ccf;        % 频率域互相关, 大小为[length(frq),length(ccr)]
        
        % 是否加窗
        switch ops.Win
            case "Hamming_1"
                taper = 0.54-0.46*cos(2*pi*ccr./max(ccr));
            case "Hamming_2"
                taper = 0.54-0.46*cos(2*pi*(ccr-min(ccr))./(max(ccr)-min(ccr)));
            case "Hamming_half"
                taper = 0.54-0.46*cos(2*pi*ccr./max(ccr));[~,tin]=max(taper);taper(1:1:tin)=1;
            case "Gauss"
                taper = exp(-1/2*((ccr-max(ccr)/2)/((max(ccr)-1)/(2*2.5))).^2);
            otherwise
                taper = 1;
        end
        ccf = ccf.*taper';

        % 检查要扫描的频段
        index= frq>=f_bound(1) & frq<=f_bound(2);
        ccf = ccf(index,:);
        frq = frq(index);

        % 如果GPU可用则将数据转移到GPU
        dispersion = zeros(length(c_range),length(ccf(:,1)));
        frq=frq(:);
        c_range=c_range(:)';   
        if canUseGPU() && strcmp(ops.GPU,"Yes")
            cb=whos("ccf",'dispersion');
            gm=gpuDevice().AvailableMemory;
            if 1*cb(1).bytes+cb(2).bytes < gm*0.8 
                frq       =gpuArray(frq);
                c_range   =gpuArray(c_range);
                ccf       =gpuArray(ccf);
                ccr       =gpuArray(ccr);
                dispersion=gpuArray(dispersion);
            end
        end

        % 开始按频率扫描
        if strcmp(ops.Fun,'J0')
            parfor (i=1:length(frq),ops.Num)
                k=2*pi*frq(i)./c_range;
                dispersion(:,i)=trapz(ccr,ccr.*ccf(i,:).'.*besselj(0,k.*ccr),1); %#ok<*PFBNS> 
            end
        elseif strcmp(ops.Fun,'J1')
            parfor (i=1:length(frq),ops.Num)
                k=2*pi*frq(i)./c_range;
                dispersion(:,i)=trapz(ccr,ccr.*ccf(i,:).'.*besselj(1,k.*ccr),1); %#ok<*PFBNS> 
            end
        elseif strcmp(ops.Fun,'H1')
            Gfs=hilbert(real(ccf));
            parfor (i=1:length(frq),ops.Num)
                k=2*pi*frq(i)./c_range;
                bjy=(besselj(0,k.*ccr)+1i*bessely(0,k.*ccr));                
                bjy(isnan(bjy))=0;bjy(isinf(bjy))=0;
                dispersion(:,i)=trapz(ccr,ccr.*Gfs(i,:).'.*bjy,1);  % ! besselh(0,2,k.*r)
            end
        elseif strcmp(ops.Fun,'H2')
            Gfs=hilbert(real(ccf));
            parfor (i=1:length(frq),ops.Num)
                k=2*pi*frq(i)./c_range;
                bjy=(besselj(0,k.*ccr)-1i*bessely(0,k.*ccr));
                bjy(isnan(bjy))=0;bjy(isinf(bjy))=0;
                dispersion(:,i)=trapz(ccr,ccr.*Gfs(i,:).'.*bjy,1);  % ! besselh(0,2,k.*r)
            end
        elseif strcmp(ops.Fun,'H2_r')
            Gfs=hilbert(real(ccf));
            parfor (i=1:length(frq),ops.Num)
                k=2*pi*frq(i)./c_range;
                bjy=(besselj(1,k.*ccr)-1i*bessely(1,k.*ccr));
                bjy(isnan(bjy))=0;bjy(isinf(bjy))=0;
                dispersion(:,i)=trapz(ccr,ccr.*Gfs(i,:).'.*bjy,1);  % ! besselh(1,2,k.*r)
            end
        elseif strcmp(ops.Fun,'FK')
            parfor (i=1:length(frq),ops.Num)
                k=2*pi*frq(i)./c_range;
                dispersion(:,i)=trapz(ccr,ccf(i,:).'.*exp(1i*k.*ccr),1);   % F-K
            end
        end

        FJ.dsp   =gather(dispersion.*sign(sum(abs(ccf),2))');  % f-c域的频散谱
        FJ.frq=gather(frq);         % 频率
        FJ.vel=gather(c_range);     % 相速度
        FJ.Fun=ops.Fun;
        FJ.Win=ops.Win;
        FJ.Rrg=gather([min(ccr) max(ccr)]);
    end


    % k滤波 @(Zhang et al., 2023)
    function CC=filtering(CC,vbound)
        arguments
            CC (1,1) struct
            vbound (1,2) {mustBeNumeric,mustBeNonnegative}=[0,0]
        end
        df=CC.freq(2)-CC.freq(1);
        if isa(CC.ccr,"cell")
            ccr = cell2mat(CC.ccr(:,1)); % 互相关对的距离
        elseif isa(CC.ccr,"double")
            ccr = CC.ccr(:,1);
        end
        for i=1:length(CC.ccf(1,:))
            if vbound(1)~=0 && 1/df/2 > ccr(i)/vbound(1)
                [pB,pA]=butter(2,ccr(i)/vbound(1)/(1/df/2),"low");
                if sum(abs(real(CC.ccf(:,i))))~=0
                    Creal=filtfilt(pB,pA,real(CC.ccf(:,i)));
                else
                    Creal=0;
                end
                if sum(abs(imag(CC.ccf(:,i))))~=0
                    Cimag=filtfilt(pB,pA,imag(CC.ccf(:,i)));
                else
                    Cimag=0;
                end
                CC.ccf(:,i)=Creal+1i*Cimag;
            end
            if vbound(2)~=0
                [pB,pA]=butter(2,ccr(i)/vbound(2)/(1/df/2),"high");
                if sum(abs(real(CC.ccf(:,i))))~=0
                    Creal=filtfilt(pB,pA,real(CC.ccf(:,i)));
                else
                    Creal=0;
                end
                if sum(abs(imag(CC.ccf(:,i))))~=0
                    Cimag=filtfilt(pB,pA,imag(CC.ccf(:,i)));
                else
                    Cimag=0;
                end
                CC.ccf(:,i)=Creal+1i*Cimag;
            end
        end
        CC.cct=fftshift( myifft_cc(CC.ccf,1+length(CC.ccf(:,1))), 1);
    end


    % 对cc.ccf按距离叠加
    function CC=stack_bin(CC,width)
        arguments
            CC (1,1) struct
            width (1,1) {mustBeNumeric,mustBeNonnegative}=0
        end
        if width==0
            return
        end
        if iscell(CC.ccr)
            R = cell2mat(CC.ccr(:,1));
        else
            R = CC.ccr(:,1);
        end
        r = R(1):width:R(end)+width/2;

        CC.ccf(CC.ccf==0)=NaN;
        ccf=zeros(length(CC.freq),length(r));
        ccr=zeros(length(r),1);
        for i=1:length(r)
            rid=R>=r(i)-width/2 & R<r(i)+width/2;
            ccf(:,i)=mean(CC.ccf(:,rid),2,"omitmissing");
            ccr(i)=mean(R(rid),1,"omitmissing");
        end
        CC.ccf=ccf;
        CC.ccr=ccr;
        CC.cct=fftshift( myifft_cc(CC.ccf,1+length(CC.ccf(:,1))), 1);

        CC.ccr(isnan(CC.ccf(1,:)),:)=[];
        CC.cct(:,isnan(CC.ccf(1,:)))=[];
        CC.ccf(:,isnan(CC.ccf(1,:)))=[];
    end

    
    function CC=inversion(FJ, CR, Mask, ops)
        arguments
            FJ (1,1) struct
            CR (:,1) {mustBeNumeric,mustBeNonnegative}
            
            Mask {mustBeNumeric}=[]
            ops.Threshold {mustBeNumeric,mustBeNonnegative}=0.01
            ops.Plt {mustBeMember(ops.Plt,["real","imag","abs"])}="abs"
            ops.GPU {mustBeMember(ops.GPU,["No","Yes"])}="No"  
            ops.Num (1,1){mustBeNumeric,mustBeNonnegative}=0;
            ops.ShowMask {mustBeMember(ops.ShowMask,["No","Yes"])}="No"  
        end
        CC=struct;

        if strcmp(ops.GPU,"Yes") && ops.Num>0
            warning('CPU + GPU parallel execution may slow down the program. It is recommended to enable only one mode.')
        end

        % FJ 是 mCCFJ.transform 函数计算得出的
        frq = FJ.frq(:);   % 频率
        vel = FJ.vel(:)';  % 相速度
        dsp = FJ.dsp;      % f-c域的频散谱, 大小为[length(frq),length(ccr)]
        fun = FJ.Fun;      % 正变换所用的基函数
        win = FJ.Win;      % 正变换所用的窗函数
        rrg = FJ.Rrg;      % 正变换的最小、最大互相关距离
        
        %
        df=frq(2)-frq(1);
        nf=length(0:df:frq(1))-1;

        % taper
        if CR(1)<rrg(1) || CR(end)>rrg(2)
            warning('Exceeding the original integration space, the result may be inaccurate!')
        end
        switch win
            case "Hamming_1"
                taper = 0.54-0.46*cos(2*pi*CR./rrg(2));
            case "Hamming_2"
                taper = 0.54-0.46*cos(2*pi*(CR-rrg(1))./(rrg(2)-rrg(1)));
            case "Hamming_half"
                taper = 0.54-0.46*cos(2*pi*CR./max(CR));[~,tin]=max(taper);taper(1:1:tin)=1;
            case "Gauss"
                taper = exp(-1/2*((CR-rrg(2)/2)/((rrg(2)-1)/(2*2.5))).^2);
            otherwise
                taper = 1;
        end
        CR=CR(:);

        % 框选感兴趣的面波频散能量团
        if isempty(Mask)
            H=figure;
            if strcmp(ops.Plt,"real")
                imagesc(frq,vel,abs(real(dsp)./max(abs(real(dsp)),[],1)));
            elseif strcmp(ops.Plt,"imag")
                imagesc(frq,vel,abs(imag(dsp)./max(abs(imag(dsp)),[],1)));
            else
                imagesc(frq,vel,abs(dsp)./max(abs(dsp),[],1));
            end
            clim([0 1]);
            xlabel('Freq. (Hz)');ylabel('Phase velocity (m/s)');
            set(gca,'YDir','normal');colormap("turbo");colorbar;
            title(gca,'Please delineate the area to be retained ...')
            roi = drawassisted(gca,'Color','r');
            title(gca,'Finished');pause(1);
            Mask = double(createMask(roi));  % 0101掩码
            Mpos  = roi.Position;            % 有值的边界线
            close(H)
        elseif Mask==1
            Mask=ones(size(dsp));
            Mpos=[];
        else
            Mask=double(Mask);
            Mpos=[];
        end


        % 根据框选的面波频散能量团提取频散曲线以及计算阵列响应函数
        if ~isempty(Mask)
            selected = (dsp .* Mask); % 注意被动源数据输入前取real，主动源取复数。嗯~也可都复数，依照问题来看
            Vr=zeros(length(frq),1);Va=zeros(length(frq),1);
            for i=1:length(frq)
                if sum(abs(selected(:,i)))==0
                    continue
                end
                % 把圈内的频散曲线给提取出来
                if strcmp(ops.Plt,"real")
                    [~,in]=max(abs(real(selected(:,i))));
                elseif strcmp(ops.Plt,"imag")
                    [~,in]=max(abs(imag(selected(:,i))));
                else
                    [~,in]=max(abs(selected(:,i)));
                end
                Vr(i)=FJ.vel(in);
                Va(i)=FJ.dsp(in,i);
            end

            % 根据频散曲线及阵列空间分布确定阵列响应函数
            ASF=zeros(length(vel),length(frq));
            parfor (i=1:length(Vr),ops.Num)
                if Vr(i)==0
                    continue;
                end
                w=2*pi*frq(i);
                k=w./vel;
                kv=w./Vr(i);

                ASF(:,i)=trapz(CR,taper.*CR.*besselj(0,kv.*CR).*besselj(0,k.*CR),1);
            end

            % 根据输入的阈值，选择ASF的有效范围
            AF1=zeros(size(ASF));
            for i=1:length(AF1(1,:))
                % 找最大值和与最大值最近的最小值的位置
                [~,in1]=max(ASF(:,i).*Mask(:,i));
                [~,in2]=findpeaks(-ASF(:,i));
                [~,in3]=min(abs(in1-in2));
                in4=in2(in3);
                din=abs(in1-in4);
                if isempty(in2)
                    AF1(:,i)=0;
                    continue
                end
                indx=(in1-3*din):1:(in1+3*din);
                indx(indx<1)=[];
                indx(indx>length(vel))=[];

                % 在区间内开始选择
                temp=smoothdata(abs(ASF(:,i)),"gaussian");
                AF1(indx,i)=temp(indx);
                Threshold=max(AF1(indx,i))*ops.Threshold;
                AF1(AF1(:,i)< Threshold,i)=0;
                AF1(AF1(:,i)>=Threshold,i)=1;
            end
            ASF_Mask=sign(AF1.*Mask.*abs(Vr'));
        else
            error("Wrong 'Mask'!")
        end

        % 绘图
        if strcmp(ops.ShowMask,"Yes")
            figure;
            if strcmp(ops.Plt,"real")
                imagesc(frq,vel,abs(real(dsp)./max(abs(real(dsp)),[],1)).*ASF_Mask);
            elseif strcmp(ops.Plt,"imag")
                imagesc(frq,vel,abs(imag(dsp)./max(abs(imag(dsp)),[],1)).*ASF_Mask);
            else
                imagesc(frq,vel,abs(dsp)./max(abs(dsp),[],1).*ASF_Mask);
            end
            hold on
            plot(frq,Vr,'w*');
            clim([0 1]);
            xlabel('Freq. (Hz)');ylabel('Phase velocity (m/s)');
            set(gca,'YDir','normal');colormap("turbo");colorbar;
            drawnow;
        end

        % 执行反变换
        if ~isempty(find(Mask==1, 1))
            selected = (dsp .* ASF_Mask); % 注意被动源数据输入前取real，主动源取复数。嗯~也可都复数，依照问题来看

            if canUseGPU() && strcmp(ops.GPU,"Yes")
                frq=gpuArray(frq);
                vel=gpuArray(vel);
                CR=gpuArray(CR);
                selected=gpuArray(selected);
            end

            INV=zeros(length(frq),length(CR));
            parfor (i=1:length(frq),ops.Num)
                if sum(abs(selected(:,i)))==0
                    continue
                end
                k=2*pi*frq(i)./vel;
                if k(1)==0
                    continue
                end

                % 反变换积分
                if strcmp(fun,'J0')
                    INV(i,:)=trapz(k,k.*selected(:,i).'.*besselj(0,k.*CR),2);
                elseif strcmp(fun,'J1')
                    INV(i,:)=trapz(k,k.*selected(:,i).'.*besselj(1,k.*CR),2);
                elseif strcmp(fun,'H1')
                    bjy=(besselj(0,k.*CR)-1i*bessely(0,k.*CR));
                    bjy(isnan(bjy))=0;bjy(isinf(bjy))=0;
                    INV(i,:)=0.5.*real(trapz(k,k.*selected(:,i).'.*bjy,2));
                elseif strcmp(fun,'H2')
                    bjy=(besselj(0,k.*CR)+1i*bessely(0,k.*CR));
                    bjy(isnan(bjy))=0;bjy(isinf(bjy))=0;
                    INV(i,:)=0.5.*real(trapz(k,k.*selected(:,i).'.*bjy,2));
                elseif strcmp(fun,'H2_r')
                    bjy=(besselj(1,k.*CR)+1i*bessely(1,k.*CR));
                    bjy(isnan(bjy))=0;bjy(isinf(bjy))=0;
                    INV(i,:)=0.5.*real(trapz(k,k.*selected(:,i).'.*bjy,2));
                elseif strcmp(fun,'FK')
                    INV(i,:)=trapz(k,k.*selected(:,i).'.*exp(-1i*k.*CR),2);
                end
            end
            INV=[zeros(nf,length(CR));INV];

            % 输出
            CC.ccf=gather(-INV./taper');
            CC.cct=fftshift( myifft_cc(CC.ccf,1+length(CC.ccf(:,1))), 1);
            CC.ccr=gather(CR);
            CC.freq=gather(0:df:max(frq))';
            CC.time=((1:length(CC.cct(:,1)))-(floor(length(CC.cct(:,1))/2) + 1))./CC.freq(end)/2;
            CC.Mask=Mask;
            CC.Mpos=Mpos;
            CC.ASF_Mask=ASF_Mask;
            CC.Vr=[zeros(nf,1);Vr];
        end
    end

    function [A,W,freq]=attenuation(CC, f_bound, r_bound, a_range, w_range, ops)
        arguments
            CC (1,1) struct
            f_bound (1,2)
            r_bound (1,2)
            a_range (1,:)
            w_range (1,:)
            ops.kind {mustBeMember(ops.kind,["J0","J1"])}="J0"
            ops.envelope {mustBeMember(ops.envelope,["No","Yes"])}="No"
            ops.ShowResults {mustBeMember(ops.ShowResults,["No","Yes"])}="No";
        end
        
        fid=CC.freq>=f_bound(1) & CC.freq<=f_bound(2);
        if iscell(CC.ccr)
            CC.ccr=cell2mat(CC.ccr(:,1));
        end
        rid=CC.ccr(:,1)>=r_bound(1) & CC.ccr(:,1)<=r_bound(2);
    
        F=CC.freq(fid);
        V=CC.Vr(fid);
        R=CC.ccr(rid,1);
        Y=real(CC.ccf(fid,rid));
        freq=F;
    
        if strcmp(ops.kind,"J1")
            v=1;
        else
            v=0;
        end
    
        A=zeros(length(F),1);W=zeros(length(F),1);
        if strcmp(ops.envelope,"No")
            parfor i=1:length(F)
                k=2*pi*F(i)/V(i);
                if isnan(k) || isinf(k); continue; end
    
                y_cal = (besselj(v,k.*R).*exp(-a_range.*R));
        
                err=zeros(length(w_range),length(a_range));
                for j=1:length(w_range)    
                    err(j,:)=sqrt(sum(abs( Y(i,:)' - w_range(j).*y_cal ).^2,1));
                end
                [~,inx,iny]=min2(err);
                A(i)=a_range(inx);
                W(i)=w_range(iny);
            end
        else
            Y1=envelope(Y')';
            parfor i=1:length(F)
                k=2*pi*F(i)/V(i);
                if isnan(k) || isinf(k); continue; end
    
                y_cal = envelope(besselj(v,k.*R).*exp(-a_range.*R));
        
                err=zeros(length(w_range),length(a_range));
                for j=1:length(w_range)    
                    err(j,:)=sqrt(sum(abs( Y1(i,:)' - w_range(j).*y_cal ).^2,1));
                end
                [~,inx,iny]=min2(err);
                A(i)=a_range(inx);
                W(i)=w_range(iny);
            end
        end
    
    
        % 交互式绘图以检查拟合效果
        if strcmp(ops.ShowResults,'Yes')
            fig = uifigure('Name', 'Interactive Plot – Check Cross‑Coherency Fitting Quality', 'Position', [500 300 600 600]);
            uibutton(fig, 'Text', 'Save','Position', [20 570 50 20],'ButtonPushedFcn', @(btn, event) saveFig(fig));
    
            ax1 = uiaxes(fig, 'Position',[50 430 500 120]);
            ax1.XLabel.String = 'Freq. (Hz)';
            ax1.YLabel.String = 'c (m/s)';
            ax1.FontSize=14;
    
            ax2 = uiaxes(fig, 'Position',[50 280 500 130]);
            ax2.XLabel.String = 'Freq. (Hz)';
            ax2.YLabel.String = '\alpha (1/m)';
            ax2.FontSize=14;
    
            ax3 = uiaxes(fig, 'Position', [50 60 500 200]);
            ax3.XLabel.String = 'Distance (m)';
            ax3.YLabel.String = 'Coherency';
            ax3.XLim = [R(1), R(end)];
            ax3.FontSize=14;
            slider = uislider(fig, ...
                'Position', [100, 40, 400, 3], ...
                'Limits', [freq(1), freq(end)], ...
                'Value', freq(1));
            slider.ValueChangedFcn = @(src, event) updatePlot(src, ax1,ax2,ax3, A,W,freq,R,V,Y);
        end
    end


    % 帮助文档
    function Help(WhichFunctionName,Language)
        arguments
            WhichFunctionName {mustBeMember(WhichFunctionName,["mCCFJ","correlate","filtering","transform","inversion","distances","stack_bin","attenuation"])}="mCCFJ"
            Language {mustBeMember(Language,["En","Zh"])}="En"
        end
        HelpFunction(WhichFunctionName,Language)
    end
end

end

%% 私有子程序
function Ut=myifft_cc(Uf,NL)
    % colum        
    Uf=gather(Uf);
    [nf,~]=size(Uf);
    Uf(isnan(Uf))=0;

    if mod(NL,2)==1 % 奇数
        Uf=[Uf;conj(Uf(nf:-1:2,:))];
    else  % 偶数
        Uf(nf,:)=real(Uf(nf,:));
        Uf=[Uf;conj(Uf(nf-1:-1:2,:))];
    end

    Ut=ifft(Uf.*nf^2,[],1,'symmetric');%
end

% 更新绘图函数
function updatePlot(src, ax1,ax2,ax3, A,W,freq,R,V,Y)
    f = src.Value;   % 获取当前滑动条值（频率）
    [~,fid]=min(abs(f-freq));
    f=freq(fid);
    src.Value=f;
    
    p1=plot(ax1,freq,V,'k-','LineWidth',1.5);hold(ax1,"on")
    plot(ax1,f,V(fid),'r.','MarkerSize',10);
    plot(ax1,f,V(fid),'ro','MarkerSize',10);hold(ax1,"off")
    legend(ax1,p1,'Extracted dispersion curve','Location','best')
    box(ax1,"on")

    p1=plot(ax2,freq,A,'k-','LineWidth',1.5);hold(ax2,"on")
    plot(ax2,f,A(fid),'r.','MarkerSize',10);
    plot(ax2,f,A(fid),'ro','MarkerSize',10);hold(ax2,"off")
    legend(ax2,p1,'Extracted attenuation curve','Location','best')
    box(ax2,"on")

    p1=plot(ax3,R,Y(fid,:),'k-','LineWidth',1.5);
    hold(ax3,"on")
    p2=plot(ax3,R,W(fid).*besselj(0,2*pi.*f./V(fid).*R).*exp(-A(fid).*R),'r--','LineWidth',1.5);
    hold(ax3,"off")
    legend(ax3,[p1 p2],'Observed','Predicted')
    box(ax3,"on")
    ax3.Title.String = ['Freq. = ', num2str(f),' Hz;  ','c = ',num2str(V(fid)),' m/s;  ','\alpha = ',num2str(A(fid)),' 1/m.'];
    ax3.Title.FontSize=16;
    ax3.Title.FontWeight="bold";
    ax3.Title.FontName='Helvetica';
end

function saveFig(fig)
    defaultName = sprintf('ShowResults_%s.fig', datestr(datetime("now"), 'yyyymmdd_HHMMSS')); %#ok<DATST>
    [file, path] = uiputfile('*.fig', 'Save as FIG files', defaultName);
    if isequal(file, 0) || isequal(path, 0)
        return; 
    end
    fullPath = fullfile(path, file);
    savefig(fig, fullPath);
end

function [value,inx,iny]=min2(A)
    % 寻找矩阵的最小值(value)及其所在行(iny)列(inx)
    [temp,in1]=min(A,[],1);
    [temp,in2]=min(temp);

    value=temp;
    inx=in2;
    iny=in1(in2);
end



function HelpFunction(WhichFunctionName,Language)
if strcmp(Language,"Zh")
    switch WhichFunctionName
        case "mCCFJ"
            disp("帮助文档 - mCCFJ")
            disp(" ")
            disp("程序版本：1.6.0")
            disp("配置要求：建议使用Matlab R2023a或更高版本")
            disp("文档用法:")
            disp("          methods(mCCFJ);")
            disp("          mCCFJ.Help('FunctionName', 'En');")
            disp("          mCCFJ.Help('FunctionName', 'Zh');")
            disp("        ")
            disp("程序作者: 杨博. 邮箱: yangb2020@mail.sustech.edu.cn")
            disp('最新版本: <a href = "https://github.com/yuanxzo/mCCFJ">见Github仓库</a>')
            disp("发布日期：2026-07-15")
            disp("All Rights Reserved.")
            disp("        ")
            disp("帮助文档 - 结尾")

        case "correlate"
            disp("帮助文档 - mCCFJ.correlate: 地震波形互相关计算程序")
            disp(" ")
            disp("用法：")
            disp("      CC=mCCFJ.correlate(wave, rloc, Fs);")
            disp("      CC=mCCFJ.correlate(wave, rloc, Fs, options);")
            disp("必要输入参数: ")
            disp("      wave    至少两道地震波形数据，大小为(npts,nsta)的数组，例如: wave=randn(npts,nsta)，其中npts为波形采样点数， nsta为台站数或波形道数")
            disp("      rloc    台站的位置，其中, rloc(1,1:nsta)为经纬度中的经度或直角坐标系里的横坐标，rloc(2,1:nsta)为经纬度中的纬度或直角坐标系里的纵坐标")
            disp("      Fs      波形数据采样率，单位Hz")
            disp("可选输入参数（options）: ")
            disp("      AX      台站的位置的类型，'latlon'表示是经纬度坐标，'xyz'表示是直角坐标，默认为'latlon'")
            disp("      CL      做相关计算时的滑动窗的长度，要求0<CL<=npts，默认CL=npts")
            disp("      OL      滑动窗的重叠长度，0<=OL<CL，默认OL=0")
            disp("      FL      对滑动窗执行傅里叶变换的长度，默认FL=CL，但是建议FL=2^nextpow2(CL)")
            disp("      NT      时间域归一化，'No'为不执行，'OneBit'为OneBit归一化，默认为'No'")
            disp("      NF      频率域归一化，'No'为不执行，'PSD'为按照该时窗内的所有台的功率谱密度的平均执行归一化，'ABS'为各台站除以自身的模的完全谱白化，默认为'No'")
            disp("      FM      所关心的最高频率，要求FM<=Fs/2，默认FM=Fs/2")
            disp("      FD      对波形的频谱做downsample的采样因子，默认为1即不做downsample。做downsample可以加快计算速度，但需要注意它对时间域互相关结果的影响")
            disp("      TP      执行傅里叶变换时的加窗类型，['No','Hann','tukeywin_5','tukeywin_10']，默认为'tukeywin_5'，即5%的tukeywin窗")
            disp("      RR      根据mCCFJ.distance函数计算得到的互相关距离信息，默认为[]，即程序会自动调用mCCFJ.distance函数计算。如果用户指定了RR，那么rloc可以为[]，AX可以不指定。")
            disp("      GPU     计算过程中是否调用GPU加快计算，'No'为不使用，'Yes'为使用，默认为'No'")
            disp("输出参数（CC，是一个结构体）: ")
            disp("      CC.acf  各台站数据的频率域自相关函数")
            disp("      CC.ccf  以站间距大小排序的频率域互相关函数")
            disp("      CC.cct  CC.ccf在时间域的结果")
            disp("      CC.ccr  第一列为互相关距离，第二列和第三列为对应互相关距离的两个台站检索，第四列为方位角，第五列和第六列为对应互相关距离的两个台站的名称")
            disp("      CC.freq 频率域互相关函数对应的频率序列")
            disp("      CC.time 时间域互相关函数对应的时间序列")
            disp("      CC.info 计算过程中保留的一些额外信息，其中，CC.info.cstack记录了各列互相关的有效叠加次数")
            disp("帮助文档 - 结尾")

        case "transform"
            disp("帮助文档 - mCCFJ.transform: 对互相关函数进行频散分析的程序")
            disp(" ")
            disp("用法：")
            disp("      FJ=mCCFJ.transform(CC, c_range, f_bound);")
            disp("      FJ=mCCFJ.transform(CC, c_range, f_bound, options);")
            disp("必要输入参数: ")
            disp("      CC      mCCFJ.correlate程序的输出")
            disp("      c_range 要成像的相速度范围，大小(1,:)，示例：c_range=linspace(200, 500, 1001); 单位m/s")
            disp("      f_bound 要成像的频率范围的上下限，大小(1,2)，示例，f_bound=[0, 80]; 单位Hz")
            disp("可选输入参数（options）: ")
            disp("      Fun     变换的类型，'J0'、'J1'、'H1'和'H2'分别为以J0、J1、H1和H2函数为基底的频率-贝塞尔变换，'FK'为F-K变换，默认为'H2'")
            disp("      Win     空间加窗的类型，['No','Hamming_1','Hamming_2','Hamming_half']，默认为'Hamming_1'")
            disp("      GPU     计算过程中是否调用GPU加快计算，'No'为不使用，'Yes'为使用，默认为'No'")
            disp("      Num     计算过程中使用线程的数量，默认为0，表示使用常规的串行计算，Num>=1时，使用parfor开启并行计算。注意，并行计算并不一定能加快计算速度。")
            disp("输出参数（FJ，是一个结构体）: ")
            disp("      FJ.frq  频率")
            disp("      FJ.vel  相速度")
            disp("      FJ.dsp  频散谱，大小为(length(FJ.vel), length(FJ.frq))的数组")
            disp("          *** 绘图时一般取频散谱的实部，但对于主动源数据一般取频散谱的绝对值")
            disp("结果展示: ")
            disp("      figure;imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(real(FJ.dsp),[],1));set(gca,'YDir','normal');clim([0 1]);");
            disp("      xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');")
            disp("帮助文档 - 结尾")

        case "inversion" 
            disp("帮助文档 - mCCFJ.inversion: 'mCCFJ.transform'的反变换程序")
            disp(" ")
            disp("用法：")
            disp("      CC=mCCFJ.inversion(FJ, RR, Mask);")
            disp("      CC=mCCFJ.inversion(FJ, RR, Mask, options);")
            disp("必要输入参数: ")
            disp("      FJ      mCCFJ.transform程序的输出。注意，对于背景噪声数据，在输入的FJ，其FJ.dsp要取实部")
            disp("      CR      距离，一般为互相关对的距离，即CC.ccr(:,1)，单位m")
            disp("      Mask    执行反变换的掩码窗，是一个与FJ.dsp相同维度大小的矩阵，全部元素由0和1组成，默认为[]，即运行程序时通过人手动画框得到")
            disp("可选输入参数（options）: ")
            disp("      Plt     频散谱绘图的类型，'real','imag'或'abs'，默认为'abs'")
            disp("      GPU     计算过程中是否调用GPU加快计算，'No'为不使用，'Yes'为使用，默认为'Yes'")
            disp("      Num     计算过程中使用线程的数量，默认为0，表示使用常规的串行计算，Num>=1时，使用parfor开启并行计算。注意，并行计算并不一定能加快计算速度。")
            disp("输出参数（CC，是一个结构体）")
            disp("帮助文档 - 结尾")

        case "attenuation"
            disp("帮助文档 - mCCFJ.attenuation: 面波衰减曲线提取程序")
            disp(" ")
            disp("用法：")
            disp("      CC=mCCFJ.attenuation(CC, f_bound, r_bound, a_range, w_range);")
            disp("      CC=mCCFJ.attenuation(CC, f_bound, r_bound, a_range, w_range, options);")
            disp("必要输入参数: ")
            disp("      CC      mCCFJ.correlate或mCCFJ.inversion程序的输出")
            disp("      f_bound 要分析的频率范围的上下限，大小(1,2)，示例，f_bound=[0, 80]; 单位Hz")
            disp("      r_bound 要分析的距离范围的上下限，大小(1,2)，示例，r_bound=[5, 50]; 单位m")
            disp("      a_range 衰减系数的网格搜索序列，大小(1,:)，示例，a_range=linspace(1e-5,1e-3,1001); 单位1/m")
            disp("      w_range 能量强度因子的网格搜索序列，大小(1,2)，示例，w_bound=linspace(0,1,1001).")
            disp("可选输入参数（options）: ")
            disp("      kind        互相干随距离的振荡类型，'J0'或'J1'，默认为'J0'，表示符合瑞雷波的振荡特征")
            disp("      envelope    是否对互相干取包络，'Yes'或'No'，默认为'No'")
            disp("      ShowResults 衰减曲线提取完成后，是否执行相干曲线拟合检查程序，'Yes'或'No'，默认为'No'")
            disp("输出参数：")
            disp("      A        衰减曲线，列向量")
            disp("      W        能量强度因子曲线，列向量")
            disp("      freq     A和W的频率序列。 ")
            disp("帮助文档 - 结尾")

        case "filtering"
            disp("帮助文档 - mCCFJ.filtering: 对互相关函数进行k滤波的程序")
            disp(" ")
            disp("用法：")
            disp("      CC=mCCFJ.filtering(CC, vbound);")
            disp("必要输入参数: ")
            disp("      CC      mCCFJ.correlate程序的输出")
            disp("      vbound  滤波的速度范围，大小(1,2)，示例：vbound=[100, 1000]; 单位m/s，意为将波速范围在小于100 m/s和大于1000 m/s的信号滤除")
            disp("输出参数（更新了CC.ccf和CC.cct的CC，是一个结构体）: ")
            disp("帮助文档 - 结尾")

        case "stack_bin"
            disp("帮助文档 - mCCFJ.stack_bin: 对'CC.ccf'按距离合并的程序")
            disp(" ")
            disp("用法：")
            disp("      CC=mCCFJ.stack_bin(CC, width);")
            disp("必要输入参数: ")
            disp("      CC      mCCFJ.correlate程序的输出")
            disp("      width   合并互相关函数的距离宽度，大小(1,1)，示例：width=10; 单位m/s，意为将互相关距离在10m范围内的所有互相关函数合并为一个")
            disp("输出参数（更新了CC.ccf和CC.cct的CC，是一个结构体）: ")
            disp("帮助文档 - 结尾")

        case "distances"
            disp("帮助文档 - mCCFJ.distances: 计算给定坐标的台站的互相关距离及方位角")
            disp(" ")
            disp("用法：")
            disp("      RR=mCCFJ.distances(rloc, AX);")
            disp("      RR=mCCFJ.distances(rloc, AX, options);")
            disp("必要输入参数: ")
            disp("      rloc    台站的位置，其中, rloc(1,:)为经纬度中的经度或直角坐标系里的横坐标，")
            disp("                               rloc(2,:)为经纬度中的纬度或直角坐标系里的纵坐标，")
            disp("                               rloc(3,:)为高程或直角坐标系里的垂向坐标（可忽略）")
            disp("      AX      台站的位置的类型，'latlon'表示是经纬度坐标，'xyz'表示是直角坐标，默认'latlon'")
            disp("可选输入参数（options）: ")
            disp("      sta_name  各台站的台站名称，默认为1:nsta")
            disp("输出参数（RR，是一个n行4列的数组）: ")
            disp("      RR(:,1) 为互相关距离，单位m；RR(:,2)为第一个台站的索引；RR(:,3)为第二个台站的索引；RR(:,4)为方位角，单位度")
            disp("帮助文档 - 结尾")

        otherwise
            disp('暂无此函数的帮助文档！')
    end
else
    switch WhichFunctionName
    case "mCCFJ"
        disp("Help Document - mCCFJ")
        disp(" ")
        disp("Program Version: 1.6.0")
        disp("Configuration Requirements: Matlab R2023a or later versions are recommended")
        disp("Document Usage:")
        disp("     methods(mCCFJ);")
        disp("     mCCFJ.Help('FunctionName', 'En');")
        disp("     mCCFJ.Help('FunctionName', 'Zh');")
        disp("        ")
        disp("Author: Bo Yang. Mail: yangb2020@mail.sustech.edu.cn")
        disp('Repository: <a href = "https://github.com/yuanxzo/mCCFJ">The mCCFJ Github repository</a>')
        disp("Release Date: 2026-07-15")
        disp("All Rights Reserved")
        disp("        ")
        disp("Help Document - End")
    case "correlate"
        disp("Help Document - mCCFJ.correlate: Seismic waveform cross-correlation calculation program")
        disp(" ")
        disp("Usage:")
        disp("      CC=mCCFJ.correlate(wave, rloc, Fs);")
        disp("      CC=mCCFJ.correlate(wave, rloc, Fs, options);")
        disp("Required Input Parameters: ")
        disp("      wave    At least two seismic waveform data, an array of size (npts,nsta), for example: wave=randn(npts,nsta), where npts is the number of waveform sampling points, nsta is the number of stations or waveform channels")
        disp("      rloc    The location of the station, where rloc(1,1:nsta) is the longitude in latitude and longitude coordinates or the x-coordinate in rectangular coordinates,")
        disp("                  rloc(2,1:nsta) is the latitude in latitude and longitude coordinates or the y-coordinate in rectangular coordinates")
        disp("                  rloc(3,1:nsta) is the elevation or vertical coordinate in rectangular coordinates (can be ignored)")
        disp("      Fs      Sampling rate of waveform data, unit Hz")
        disp("Optional Input Parameters (options): ")
        disp("      AX      The type of station location, 'latlon' means it is latitude and longitude coordinates, 'xyz' means it is rectangular coordinates, default 'latlon'")
        disp("      CL      The length of the sliding window for correlation calculation, requires 0<CL<=npts, default CL=npts")
        disp("      OL      The overlap length of the sliding window, 0<=OL<CL, default OL=0")
        disp("      FL      The length for performing Fourier transform on sliding windows, default FL=CL. A recommended effective value is FL=2^nextpow2(CL) to accelerate computation, or FL=2*CL-1 to align with the 'xcorr' results.")
        disp("      NT      Time domain normalization, 'No' means not executed, 'OneBit' means OneBit normalization, default 'No'")
        disp("      NF      Frequency domain normalization, 'No' means not executed, 'PSD' means normalized by the average power spectral density of all stations in this time window, 'ABS', complete spectral whitening, means each station divided by its own modulus, default 'No'")
        disp("      FM      The maximum frequency of interest, requires FM<=Fs/2, default FM=Fs/2")
        disp("      FD      Downsample factor for the waveform spectrum, default 1 means no downsample. Downsampling can speed up the calculation, but it is necessary to pay attention to its impact on the time domain cross-correlation result")
        disp("      TP      The taper of Fourier transformation, ['No','Hann','tukeywin_5','tukeywin_10'], default 'tukeywin_5', i.e., 5% tukeywin window")
        disp("      RR      The inter-station distance information calculated by the 'mCCFJ.distance' function, default [] means the program will automatically call the 'mCCFJ.distance' function to calculate. If 'RR' is specified, 'rloc' can be [] and 'AX' can be the default.")
        disp("      GPU     Whether to use GPU to speed up the calculation during the calculation process, 'No' means not using, 'Yes' means using, default 'No'")
        disp("Output Parameters (CC, is a structure): ")
        disp("      CC.acf  Frequency domain autocorrelation function of each station data")
        disp("      CC.ccf  Frequency domain cross-correlation function sorted by inter-station distance")
        disp("      CC.cct  Time-domain results of CC.ccf")
        disp("      CC.ccr  The 1st column is the inter-station distance, the 2nd and 3rd columns are the station retrieval, and the 4th column is the azimuth, the 5th and 6th columns are the station names")
        disp("      CC.freq Frequency sequence corresponding to the frequency domain cross-correlation function")
        disp("      CC.time Time sequence corresponding to the time domain cross-correlation function")
        disp("      CC.info Some additional information reserved during the calculation process, where CC.info.cstack records the effective stacking times of each column of CC.ccf")
        disp("Help Document - End")

    case "transform"
        disp("Help Document - mCCFJ.transform: Dispersion analysis program for cross-correlation function")
        disp(" ")
        disp("Usage:")
        disp("      FJ=mCCFJ.transform(CC, c_range, f_bound);")
        disp("      FJ=mCCFJ.transform(CC, c_range, f_bound, options);")
        disp("Required Input Parameters: ")
        disp("      CC      Output of 'mCCFJ.correlate' program")
        disp("      c_range Range of phase velocity to be imaged, size (1,:), example: c_range=linspace(200, 500, 1001); unit m/s")
        disp("      f_bound Upper and lower limits of the frequency range to be imaged, size (1,2), example: f_bound=[0, 80]; unit Hz")
        disp("Optional Input Parameters (options): ")
        disp("      Fun     Type of transformation, 'J0', 'J1', 'H1' and 'H2' are frequency-Bessel transformations with J0, J1, H1 and H2 functions as bases respectively, 'FK' is F-K transformation, default 'H2'")
        disp("      Win     Type of spatial windowing, ['No','Hamming_1','Hamming_2','Hamming_half'], default 'Hamming_1'")
        disp("      GPU     Whether to use GPU to speed up the calculation during the calculation process, 'No' means not using, 'Yes' means using, default 'No'")
        disp("      Num     The number of threads used in the computation process. The default is 0, indicating conventional serial computation. When Num >= 1, parallel computation is enabled using 'parfor'. Note that 'parfor' does not necessarily speed up the computation.")
        disp("Output Parameters (FJ, is a structure): ")
        disp("      FJ.frq  Frequency")
        disp("      FJ.vel  Phase velocity")
        disp("      FJ.dsp  Dispersion spectrum, size (length(FJ.vel), length(FJ.frq)) array")
        disp("          *** When plotting, the real part of the dispersion spectrum is generally taken, but for active source data, the absolute value of the dispersion spectrum is generally taken. (Non absolute)")
        disp("Result Display: ")
        disp("      figure;imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(real(FJ.dsp),[],1));set(gca,'YDir','normal');clim([0 1]);");
        disp("      xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');")
        disp("Help Document - End")

    case "inversion"
        disp("Help Document - mCCFJ.inversion: Inverse transformation program for 'mCCFJ.transform'")
        disp(" ")
        disp("Usage:")
        disp("      CC=mCCFJ.inversion(FJ, RR, Mask);")
        disp("      CC=mCCFJ.inversion(FJ, RR, Mask, options);")
        disp("Required Input Parameters: ")
        disp("      FJ     Output of 'mCCFJ.transform' program. Note: For ambient noise data, the real part of FJ.dsp should be taken as input")
        disp("      CR     Distances, generally the inter-station distances from cross-correlation pairs (i.e., CC.ccr(:,1)), in meters")
        disp("      Mask   Mask window for performing the inverse transformation. It is a matrix of the same dimensions as FJ.dsp, with all elements consisting of 0 and 1. Default is [], meaning it will be manually drawn during program execution")
        disp("Optional Input Parameters (options): ")
        disp("      Plt    Type of dispersion spectrum plot: 'real', 'imag', or 'abs'. Default is 'abs'")
        disp("      GPU    Whether to use GPU acceleration during computation: 'No' for disable, 'Yes' for enable. Default is 'Yes'")
        disp("      Num     The number of threads used in the computation process. The default is 0, indicating conventional serial computation. When Num >= 1, parallel computation is enabled using 'parfor'. Note that 'parfor' does not necessarily speed up the computation.")
        disp("Output Parameter (CC, a structure)")
        disp("Help Document - End")

    case "attenuation"
        disp("Help document - mCCFJ.attenuation: Surface-wave attenuation curve extraction program")
        disp(" ")
        disp("Usage:")
        disp("      CC=mCCFJ.attenuation(CC, f_bound, r_bound, a_range, w_range);")
        disp("      CC=mCCFJ.attenuation(CC, f_bound, r_bound, a_range, w_range, options);")
        disp("Required Input Parameters: ")
        disp("      CC      Output of 'mCCFJ.correlate' or 'mCCFJ.inversion' program")
        disp("      f_bound Lower and upper limits of the frequency range to analyze, size (1,2), e.g., f_bound=[0, 80]; unit Hz")
        disp("      r_bound Lower and upper limits of the distance range to analyze, size (1,2), e.g., r_bound=[5, 50]; unit m")
        disp("      a_range Grid search sequence for attenuation coefficient, size (1,:), e.g., a_range=linspace(1e-5,1e-3,1001); unit 1/m")
        disp("      w_range Grid search sequence for energy intensity factor, size (1,2), e.g., w_bound=linspace(0,1,1001).")
        disp("Optional Input Parameters (options): ")
        disp("      kind        Oscillation type of the coherency with distance, 'J0' or 'J1', default 'J0', representing Rayleigh-wave oscillation characteristics.")
        disp("      envelope    Whether to take the envelope for the coherency, 'Yes' or 'No', default 'No'")
        disp("      ShowResults Whether to run the coherency curve fitting check after attenuation extraction, 'Yes' or 'No', default 'No'")
        disp("Output parameters:")
        disp("      A        Attenuation curve, column vector, unit 1/m")
        disp("      W        Energy intensity factor curve, column vector")
        disp("      freq     Frequency vector corresponding to A and W, unit Hz.")
        disp("Help Document - End")

    case "filtering"
        disp("Help Document - mCCFJ.filtering: k-filtering program for cross-correlation function")
        disp(" ")
        disp("Usage:")
        disp("      CC=mCCFJ.filtering(CC, vbound);")
        disp("Required Input Parameters: ")
        disp("      CC      Output of 'mCCFJ.correlate' program")
        disp("      vbound  Filtering velocity range, size (1,2), example: vbound=[100, 1000]; unit m/s, meaning to filter out signals with wave speeds less than 100 m/s and greater than 1000 m/s")
        disp("Output Parameters (updated CC.ccf and CC.cct, CC is a structure): ")
        disp("Help Document - End")

    case "stack_bin"
        disp("Help Document - mCCFJ.stack_bin: Program for Merging 'CC.ccf' by Distance")
        disp(" ")
        disp("Usage:")
        disp("      CC=mCCFJ.stack_bin(CC, width);")
        disp("Required input parameters:")
        disp("      CC      Output of 'mCCFJ.correlate' program")
        disp("      width   The width of distance bin for stacking correlations, size (1,1), example: width=10; unit m, which means merging all cross-correlation functions within a range of 10 m into one")
        disp("Output Parameters (updated CC.ccf and CC.cct, CC is a structure): ")
        disp("Help Document - End")

    case "distances"
        disp("Help Document - mCCFJ.distances: Calculate the inter-station distance and azimuth of the given coordinates")
        disp(" ")
        disp("Usage:")
        disp("      RR=mCCFJ.distances(rloc, AX);")
        disp("      RR=mCCFJ.distances(rloc, AX, options);")
        disp("Required Input Parameters: ")
        disp("      rloc    The location of the station, where rloc(1,:) is the longitude in latitude and longitude coordinates or the x-coordinate in rectangular coordinates,")
        disp("                  rloc(2,:) is the latitude in latitude and longitude coordinates or the y-coordinate in rectangular coordinates,")
        disp("                  rloc(3,:) is the elevation or vertical coordinate in rectangular coordinates (can be ignored)")
        disp("      AX      The type of station location, 'latlon' means it is latitude and longitude coordinates, 'xyz' means it is rectangular coordinates, default AX=0")
        disp("Optional Input Parameters (options): ")
        disp("      sta_name  The station names of each station, default 1: nsta")
        disp("Output Parameters (RR is a n by 4 array): ")
        disp("      RR(:,1) is the inter-station distance, unit m; RR(:,2) is the index of the first station; RR(:,3) is the index of the second station; RR(:,4) is the azimuth, unit degree")
        disp("Help Document - End")

    otherwise
        disp('No help document for this function yet!')
    end
end
end
