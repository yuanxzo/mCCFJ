classdef mCCFJ
% A MATLAB Package for calculating seismic ambient noise cross-correlation and frequency-Bessel transformation.
% 
% In this package, these are two main function, mCCFJ.correlate and mCCFJ.transform.
%    mCCFJ.correlate is used to calculate the cross-correlation or
%       cross-coherency function of seismic waveforms. To ensure the
%       efficiency of calculation, the function is calculated in the
%       frequency domain, and the cross-correlation between any two
%       stations in the same time window is calculated in the way of matrix
%       parallelism. If necessary, GPU acceleration can be used.
%    mCCFJ.transform is used for dispersion analysis of cross-correlation
%       function. We provide a variety of frequency-wavenumber domain
%       transformation methods to deal with different data, which is up to
%       you. In order to ensure the efficiency of computing, GPU
%       acceleration can be used when necessary.
% To_begin_with = 'mCCFJ.Help'

% Accessible program: called by using mCCFJ.functionname
methods (Static)

    % Correlate seismic waveforms
    function CC=correlate(wave,rloc,Fs,ops)
        arguments
            wave (:,:) {mustBeNumeric}
            rloc (:,:) {mustBeNumeric} = [1:1:length(wave(1,:));zeros(1,length(wave(1,:)));zeros(1,length(wave(1,:)))]
            Fs (1,1) {mustBeNumeric,mustBePositive} =1
            ops.AX {mustBeMember(ops.AX,["latlon","xyz"])}="latlon"
            ops.CL (1,1) {mustBeInteger,mustBePositive} = length(wave(:,1))
            ops.OL (1,1) {mustBeInteger,mustBeNonnegative} = 0
            ops.NT {mustBeMember(ops.NT,["No","OneBit"])}="No"
            ops.NF {mustBeMember(ops.NF,["No","PSD","ABS"])}="No"
            ops.FM (1,1) {mustBeNumeric,mustBePositive} = Fs/2
            ops.FD {mustBeInteger,mustBePositive}=1    % Frequency point dilution
            ops.TP {mustBeMember(ops.TP,["No","Hann","tukeywin_5","tukeywin_10"])}="tukeywin_5"
            ops.RR (:,4) {mustBeNumeric} = []
            ops.GPU{mustBeMember(ops.GPU,["No","Yes"])}="No"
        end
        CC=struct;
    
        % 可选参数解码
        [npts,nsta]=size(wave,[1 2]);
        if nsta<=1
            error("At least two seismic waveform data are required!")
        end
        if npts<=1
            error("The number of sampling points in the waveform data is less than 1!")
        end
        if ops.CL > npts
            error("The length of the sliding window for correlation calculation exceeds the total length of the waveform, ending calculation!");
        end
        if ops.OL > ops.CL
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
            warning("The maximum frequency of interest is greater than or equal to half the sampling rate, setting FM=Fs/2");
            ops.FM = Fs/2;
        end
    
        % 确定互相关滑动窗的起始点的位置及个数
        index = 1 : ops.CL-ops.OL : npts-ops.CL+1;
        if npts - index(end) + 1 < ops.CL
            index(end)=[];
        end
        num_of_win = length(index);
        len_of_win = ops.CL;
    
        % 计算所有可能的台站对的站间距离
        if isempty(ops.RR)
            ccr=mCCFJ.distances(rloc,ops.AX); % ccr(:,1) 为距离, ccr(:,2) 为第一个台站的索引, cc(:,3) 为第二个台站的索引, cc(:,4) 为方位角
        else
            ccr=ops.RR;
        end
        num_of_ccr = length(ccr(:,1));
    
        % One-bit
        if strcmp(ops.NT,"OneBit")
            wave = sign(wave);
        end
    
        % 频率序列
        freq = Fs/(ops.CL)*(0:1:ops.CL/2);
        finx = freq<=ops.FM;
        freq = freq(finx);
        freq = freq(1:ops.FD:end);
        
        % 分配内存
        acf = zeros(length(freq),nsta);
        ccf = zeros(length(freq),num_of_ccr);
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
            otherwise
                error("Unknown taper type, ending calculation!");
        end

        % 检查GPU是否可用
        if canUseGPU() && strcmp(ops.GPU,"Yes")
            cb=whos("ccf",'wave');
            gm=gpuDevice().AvailableMemory;
            if 2*cb(1).bytes+cb(2).bytes < gm*0.8 
                wave=gpuArray(wave);
                acf=gpuArray(acf);
                ccf=gpuArray(ccf);
                ccr=gpuArray(ccr);
                index=gpuArray(index); 
                len_of_win=gpuArray(len_of_win);
                num_of_win=gpuArray(num_of_win);
                taper=gpuArray(taper);
                finx=gpuArray(finx);
                astack=gpuArray(astack);
                cstack=gpuArray(cstack);

                GPUtype=0;
            elseif 2*cb(1).bytes < gm*0.8 
                acf=gpuArray(acf);
                ccf=gpuArray(ccf);
                ccr=gpuArray(ccr);
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
                disp("Attempted to use the GPU, but graphics memory is insufficient. CPU is running.")
            end
        else
            GPUtype=0;
        end
    

        % 逐个滑动窗计算互相关, 平均为最终结果
        for i=1:num_of_win
            
            wave_seg = wave(index(i):index(i)+len_of_win-1,:);

            % 如果有些道的数据含有0过多则不用此道
            temp = sum(sign(abs(wave_seg)),1);
            inx0 = temp<=len_of_win-5;
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
            wave_seg = fft(wave_seg.*taper,[],1)/len_of_win;
            wave_seg = wave_seg(1:len_of_win/2+1,:);
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

            % 互相关 和 谱白化
            st1 = int64(ccr(:,2)); st2 = int64(ccr(:,3));
            cc_temp = conj(wave_seg(:,st1)).*wave_seg(:,st2);
            switch ops.NF
                case "No"   % 直接互相关
                    cc_norm = 1;
                case "PSD"  % PSD归一化
                    cc_norm = mean(temp(:,inx1),2); % 当前时窗的PSD
                case "ABS"  % 绝对值归一化，完全的谱白化
                    cc_norm = sqrt(temp(:,st1)).*sqrt(temp(:,st2));
                otherwise
                    error("Unknown normalization type, ending calculation!");
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

        % 不符合要求的强行置为0
        inx4 = astack==0; acf(:,inx4) = 0;
        inx4 = cstack==0; ccf(:,inx4) = 0;
        
        % 计算时间域互相关
        cct = fftshift( myifft(ccf), 1);

        % 输出结果    
        CC.acf = gather(acf);
        CC.ccf = gather(ccf);
        CC.cct = gather(cct);
        CC.ccr = gather(ccr);
        CC.time = gather(linspace(-ops.CL/Fs/2/ops.FD,ops.CL/Fs/2/ops.FD,length(cct(:,1))));
        CC.freq = gather(freq');
        CC.info.astack = gather(astack);
        CC.info.cstack = gather(cstack);
        CC.info.num_of_win = gather(num_of_win);
    end

    % 计算互相关距离
    function RR=distances(rloc,AX)
        arguments
            rloc (:,:) {mustBeNumeric}
            AX   {mustBeMember(AX,["latlon","xyz"])}
        end
        if length(rloc(:,1))==2
            rloc(3,:)=0;
        end
        nsta=length(rloc(1,:));
        RR = zeros((nsta*(nsta-1)/2),4);  % RR(:,1) 为距离, RR(:,2) 为第一个台站的索引, RR(:,3) 为第二个台站的索引, RR(:,4) 为方位角
        nccr = 0;
        if strcmp(AX,"xyz")     % 直角坐标
            for i = 1:nsta-1
                s1=rloc(:,i);
                s2=rloc(:,i+1:nsta);
                RR(nccr+1:nccr+nsta-i,1)=sqrt(sum((s1-s2).^2,1));
                RR(nccr+1:nccr+nsta-i,2)=i;
                RR(nccr+1:nccr+nsta-i,3)=i+1:nsta;
                RR(nccr+1:nccr+nsta-i,4)=mod(450-rad2deg(atan2(s2(2,:)-s1(2), s2(1,:)-s1(1))), 360);

                nccr=nccr+nsta-i;
            end
        elseif strcmp(AX,"latlon") % 经纬度坐标
            for i = 1:nsta-1
                s1=rloc(:,i);
                for j = i+1:nsta
                    s2=rloc(:,j);
                    nccr=nccr+1;
                    [dist,  az]=distance(s1(2), s1(1), s2(2), s2(1), referenceEllipsoid('WGS 84'));
                    RR(nccr,1)=sqrt(dist.^2+(s2(3)-s1(3)).^2);
                    RR(nccr,2)=i;
                    RR(nccr,3)=j;
                    RR(nccr,4)=az;      % 前一个台作观察后一个台的方位角
                end
            end
        end
    
        [~, order] = sort(RR(:,1));
        RR=RR(order,:);
    end


    % 频率-贝塞尔变换
    function FJ=transform(CC, c_range, f_bound, ops)
        arguments
            CC (1,1) struct
            c_range {mustBeNumeric,mustBeNonzero}
            f_bound {mustBeNumeric,mustBeNonnegative}
            ops.Fun {mustBeMember(ops.Fun,["J0","J1","H1","H2","FK"])}="H2"
            ops.Win {mustBeMember(ops.Win,["No","Hamming_1","Hamming_2","Hamming_half"])}="Hamming_1"
            ops.GPU {mustBeMember(ops.GPU,["No","Yes"])}="Yes"            
        end
        FJ=struct;

        % CC 是 correlate 函数计算得出的, 如果不是，按下面形式合成CC也可
        frq = CC.freq;     % 频率域互相关依赖的频率序列
        ccr = CC.ccr(:,1); % 互相关对的距离
        ccf = CC.ccf;      % 频率域互相关, 大小为[length(frq),length(ccr)]
        
        % 是否加窗
        if strcmp(ops.Win,"Hamming_1")
            taper = 0.54-0.46*cos(2*pi*ccr./max(ccr));
        elseif strcmp(ops.Win,"Hamming_2")
            taper = 0.54-0.46*cos(2*pi*(ccr-min(ccr))./(max(ccr)-min(ccr)));
        elseif strcmp(ops.Win,"Hamming_half")
            taper = 0.54-0.46*cos(2*pi*ccr./max(ccr));
            taper(1:1:length(ccr)/2)=1;
        else
            taper = 1;
        end
        ccf = ccf.*taper(:)';

        % 检查要扫描的频段
        index= frq>=f_bound(1) & frq<=f_bound(2);
        ccf = ccf(index,:);
        frq = frq(index);

        % 如果GPU可用则将数据转移到GPU
        dispersion = zeros(length(c_range),length(ccf(:,1)));
        c_range=c_range(:)';
        frq=frq(:);
        ccr=ccr(:);
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
            for i=1:length(frq)
                k=2*pi*frq(i)./c_range;
                dispersion(:,i)=trapz(ccr,ccr.*ccf(i,:).'.*besselj(0,k.*ccr),1); %#ok<*PFBNS> 
            end
        elseif strcmp(ops.Fun,'J1')
            for i=1:length(frq)
                k=2*pi*frq(i)./c_range;
                dispersion(:,i)=trapz(ccr,ccr.*ccf(i,:).'.*besselj(1,k.*ccr),1); %#ok<*PFBNS> 
            end
        elseif strcmp(ops.Fun,'H1')
            Gfs=hilbert(real(ccf));
            for i=1:length(frq)
                k=2*pi*frq(i)./c_range;
                bjy=(besselj(0,k.*ccr)+1i*bessely(0,k.*ccr));
                if isinf(bjy(1,1)) || isnan(bjy(1,1))
                    bjy(1,:)=0;
                end
                dispersion(:,i)=trapz(ccr,ccr.*Gfs(i,:).'.*bjy,1);  % ! besselh(0,2,k.*r)
            end
        elseif strcmp(ops.Fun,'H2')
            Gfs=hilbert(real(ccf));
            for i=1:length(frq)
                k=2*pi*frq(i)./c_range;
                bjy=(besselj(0,k.*ccr)-1i*bessely(0,k.*ccr));
                if isinf(bjy(1,1)) || isnan(bjy(1,1))
                    bjy(1,:)=0;
                end
                dispersion(:,i)=trapz(ccr,ccr.*Gfs(i,:).'.*bjy,1);  % ! besselh(0,2,k.*r)
            end
        elseif strcmp(ops.Fun,'FK')
            for i=1:length(frq)
                k=2*pi*frq(i)./c_range;
                dispersion(:,i)=trapz(ccr,ccf(i,:).'.*exp(1i*k.*ccr),1);   % F-K
            end
        end

        FJ.dsp=gather(dispersion.*sign(sum(abs(ccf),2))');  % 频散谱
        FJ.frq=gather(frq);         % 频率
        FJ.vel=gather(c_range);     % 相速度
        FJ.Fun=ops.Fun;
        FJ.Win=ops.Win;
        FJ.Rrg=[min(ccr) max(ccr)];
    end


    % k滤波
    function CC=filtering(CC,vbound)
        arguments
            CC (1,1) struct
            vbound (1,2) {mustBeNumeric,mustBeNonnegative}=[0,0]
        end
        df=CC.freq(2)-CC.freq(1);
        for i=1:length(CC.ccf(1,:))
            if vbound(1)~=0 && 1/df/2 > CC.ccr(i,1)/vbound(1)
                [pB,pA]=butter(2,CC.ccr(i,1)/vbound(1)/(1/df/2),"low");
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
                [pB,pA]=butter(2,CC.ccr(i,1)/vbound(2)/(1/df/2),"high");
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
        CC.cct=fftshift( myifft(CC.ccf), 1);
    end


    % 帮助文档
    function Help(WhichFunctionName,Language)
        arguments
            WhichFunctionName {mustBeMember(WhichFunctionName,["mCCFJ","correlate","filtering","transform","distances"])}="mCCFJ"
            Language {mustBeMember(Language,["En","Zh"])}="En"
        end
        HelpFunction(WhichFunctionName,Language)
    end
end

end
%% 私有子程序

function Ut=myifft(Uf)
    % colum    
    Uf=gather(Uf);
    [nf,~]=size(Uf);

    Uf=Uf.*nf;
    Uf(nf,:)=real(Uf(nf,:));
    Uf=[Uf;conj(Uf(nf-1:-1:2,:))];

    Ut=ifft(Uf,[],1,'symmetric');%
end


function HelpFunction(WhichFunctionName,Language)
if strcmp(Language,"Zh")
    switch WhichFunctionName
        case "mCCFJ"
            disp("帮助文档 - mCCFJ")
            disp(" ")
            disp("程序版本：1.2.0")
            disp("配置要求：建议使用Matlab R2022b或更新版本")
            disp("文档用法:")
            disp("          methods(mCCFJ);")
            disp("          mCCFJ.Help('FunctionName', 'En');")
            disp("          mCCFJ.Help('FunctionName', 'Zh');")
            disp("        ")
            disp("程序作者: 杨博. Mail: yangb2020@mail.sustech.edu.cn")
            disp('最新版本: <a href = "https://github.com/yuanxzo/mCCFJ">见Github仓库</a>')
            disp("发布日期：2025-05-20")
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
            disp("      NT      时间域归一化，'No'为不执行，'OneBit'为OneBit归一化，默认为'No'")
            disp("      NF      频率域归一化，'No'为不执行，'PSD'为按照该时窗内的所有台的功率谱密度的平均执行归一化，'ABS'为各台站除以自身的模的完全谱白化，默认为'No'")
            disp("      FM      所关心的最高频率，要求FM<=Fs/2，默认FM<=Fs/2")
            disp("      FD      对波形的频谱做downsample的采样因子，默认为1即不做downsample。做downsample可以加快计算速度，但需要注意它对时间域互相关结果的影响")
            disp("      TP      执行傅里叶变换时的加窗类型，['No','Hann','tukeywin_5','tukeywin_10']，默认为'tukeywin_5'，即5%的tukeywin窗")
            disp("      RR      根据mCCFJ.distance函数计算得到的互相关距离信息，默认为[]，即程序会自动调用mCCFJ.distance函数计算")
            disp("      GPU     计算过程中是否调用GPU加快计算，'No'为不使用，'Yes'为使用，默认为'No'")
            disp("输出参数（CC，是一个结构体）: ")
            disp("      CC.acf  各台站数据的频率域自相关函数")
            disp("      CC.ccf  以站间距大小排序的频率域互相关函数")
            disp("      CC.cct  以站间距大小排序的时间域互相关函数")
            disp("      CC.ccr  第一列为互相关距离，第二列和第三列为对应互相关距离的台站检索，第四列为方位角")
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
            disp("      GPU     计算过程中是否调用GPU加快计算，'No'为不使用，'Yes'为使用，默认为'Yes'")
            disp("输出参数（FJ，是一个结构体）: ")
            disp("      FJ.frq  频率")
            disp("      FJ.vel  相速度")
            disp("      FJ.dsp  频散谱，大小为(length(FJ.vel), length(FJ.frq))的数组")
            disp("          *** 绘图时一般取频散谱的实部，但对于主动源数据一般取频散谱的绝对值")
            disp("结果展示: ")
            disp("      figure;imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(abs(real(FJ.dsp)),[],1));set(gca,'YDir','normal');clim([0 1]);");
            disp("      xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');")
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

        case "distances"
            disp("帮助文档 - mCCFJ.distances: 计算给定坐标的台站的互相关距离及方位角")
            disp(" ")
            disp("用法：")
            disp("      RR=mCCFJ.distances(rloc,AX);")
            disp("必要输入参数: ")
            disp("      rloc    台站的位置，其中, rloc(1,:)为经纬度中的经度或直角坐标系里的横坐标，")
            disp("                               rloc(2,:)为经纬度中的纬度或直角坐标系里的纵坐标，")
            disp("                               rloc(3,:)为高程或直角坐标系里的垂向坐标（可忽略）")
            disp("      AX      台站的位置的类型，'latlon'表示是经纬度坐标，'xyz'表示是直角坐标，默认'latlon'")
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
        disp("Program Version: 1.2.0")
        disp("Configuration Requirements: Matlab R2022b or later versions are recommended")
        disp("Document Usage:")
        disp("     methods(mCCFJ);")
        disp("     mCCFJ.Help('FunctionName', 'En');")
        disp("     mCCFJ.Help('FunctionName', 'Zh');")
        disp("        ")
        disp("Author: Bo Yang. Mail: yangb2020@mail.sustech.edu.cn")
        disp('Repository: <a href = "https://github.com/yuanxzo/mCCFJ">The mCCFJ Github repository</a>')
        disp("Release Date: 2025-05-20")
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
        disp("      NT      Time domain normalization, 'No' means not executed, 'OneBit' means OneBit normalization, default 'No'")
        disp("      NF      Frequency domain normalization, 'No' means not executed, 'PSD' means normalized by the average power spectral density of all stations in this time window, 'ABS', complete spectral whitening, means each station divided by its own modulus, default 'No'")
        disp("      FM      The maximum frequency of interest, requires FM<=Fs/2, default FM<=Fs/2")
        disp("      FD      Downsample factor for the waveform spectrum, default 1 means no downsample. Downsampling can speed up the calculation, but it is necessary to pay attention to its impact on the time domain cross-correlation result")
        disp("      TP      The taper of Fourier transformation, ['No','Hann','tukeywin_5','tukeywin_10'], default 'tukeywin_5', i.e., 5% tukeywin window")
        disp("      RR      The inter-station distance information calculated by the 'mCCFJ.distance' function, default [] means the program will automatically call the 'mCCFJ.distance' function to calculate")
        disp("      GPU     Whether to use GPU to speed up the calculation during the calculation process, 'No' means not using, 'Yes' means using, default 'No'")
        disp("Output Parameters (CC, is a structure): ")
        disp("      CC.acf  Frequency domain autocorrelation function of each station data")
        disp("      CC.ccf  Frequency domain cross-correlation function sorted by inter-station distance")
        disp("      CC.cct  Time domain cross-correlation function sorted by inter-station distance")
        disp("      CC.ccr  The first column is the inter-station distance, the second and third columns are the station retrieval corresponding to the inter-station distance, and the fourth column is the azimuth")
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
        disp("      GPU     Whether to use GPU to speed up the calculation during the calculation process, 'No' means not using, 'Yes' means using, default 'Yes'")
        disp("Output Parameters (FJ, is a structure): ")
        disp("      FJ.frq  Frequency")
        disp("      FJ.vel  Phase velocity")
        disp("      FJ.dsp  Dispersion spectrum, size (length(FJ.vel), length(FJ.frq)) array")
        disp("          *** When plotting, the real part of the dispersion spectrum is generally taken, but for active source data, the absolute value of the dispersion spectrum is generally taken. (Non absolute)")
        disp("Result Display: ")
        disp("      figure;imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(abs(real(FJ.dsp)),[],1));set(gca,'YDir','normal');clim([0 1]);");
        disp("      xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');")
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
        
    case "distances"
        disp("Help Document - mCCFJ.distances: Calculate the inter-station distance and azimuth of the given coordinates")
        disp(" ")
        disp("Usage:")
        disp("      RR=mCCFJ.distances(rloc, AX);")
        disp("Required Input Parameters: ")
        disp("      rloc    The location of the station, where rloc(1,:) is the longitude in latitude and longitude coordinates or the x-coordinate in rectangular coordinates,")
        disp("                  rloc(2,:) is the latitude in latitude and longitude coordinates or the y-coordinate in rectangular coordinates,")
        disp("                  rloc(3,:) is the elevation or vertical coordinate in rectangular coordinates (can be ignored)")
        disp("      AX      The type of station location, 'latlon' means it is latitude and longitude coordinates, 'xyz' means it is rectangular coordinates, default AX=0")
        disp("Output Parameters (RR is a n by 4 array): ")
        disp("      RR(:,1) is the inter-station distance, unit m; RR(:,2) is the index of the first station; RR(:,3) is the index of the second station; RR(:,4) is the azimuth, unit degree")
        disp("Help Document - End")

    otherwise
        disp('No help document for this function yet!')
    end
end
end
