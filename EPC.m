function[M,m,phaseCongruency,or]=EPC(im)

    minWaveLength = 1;
    sigmaOnf = 5;
    k = 0.5;

    epsilon=.0001;
    imagefft=single(fft2(im));
    sze=size(imagefft);
    rows=sze(1);
    cols=sze(2);
    zero=single(zeros(sze));
    covx2=zero;
    covy2=zero;
    covxy=zero;


    x=single(ones(rows,1)*(-cols/2:(cols/2-1))/(cols/2));
    y=single((-rows/2:(rows/2-1))'*ones(1,cols)/(rows/2));
    radius=sqrt(x.^2+y.^2);
    radius(round(rows/2+1),round(cols/2+1))=1;   
    clear x;clear y;
    
    Energy_ThisOrient=zero;
    

    [X,Y] = meshgrid(-minWaveLength:minWaveLength,-minWaveLength:minWaveLength);   
    LoG = -1/(sigmaOnf)^2 / pi * exp((-(X.^2+Y.^2))/(2*sigmaOnf^2)).*(1-((X.^2+Y.^2))/(2*(sigmaOnf)^2));
    h1 = [- 1, 0, 1; - 2, 0, 2; - 1, 0, 1];
    h2 = [- 1, - 2, - 1; 0, 0, 0; 1, 2, 1];


%% 奇对称 虚部
   GX = imfilter(im,h1,'replicate','same');
   GY = imfilter(im,h2,'replicate','same');      
   GM = sqrt(GX.^2 + GY.^2); 
   %figure;imshow(mapminmax(GM,0,1));
   
%% 偶对称 实部   
    Ln = imfilter(im,LoG,'replicate','same');
    %figure;imshow(mapminmax(Ln,0,1));
    
    filter = ((-((radius)).^2)/(2*(0.55)^2));
    filter(round(rows/2+1),round(cols/2+1))=0;
    filter=fftshift(filter);
    ifftFilt=real(ifft2(filter))*sqrt(rows*cols);
    ifftFilterArray=single(ifftFilt);
    clear ifftFilt; clear GX; clear GY; clear h1; clear h2; clear LoG; clear X; clear Y;
    
    EO = complex(Ln, GM);
    
    clear EOfft;
    An=abs(EO);
    sumAn_ThisOrient=single(An);
    sumE_ThisOrient=single(real(EO));
    sumO_ThisOrient=single(imag(EO));

    EM_n=sum(sum(filter.^2));
    clear filter;
    
    XEnergy=sqrt(sumE_ThisOrient.^2+sumO_ThisOrient.^2)+epsilon;
    MeanE=sumE_ThisOrient./XEnergy;
    MeanO=sumO_ThisOrient./XEnergy;
    
    E=real(EO);O=imag(EO);
    Energy_ThisOrient=Energy_ThisOrient+E.*MeanE+O.*MeanO-abs(E.*MeanO-O.*MeanE);

    
    medianE2n=median(reshape(An.^2,1,rows*cols));
    meanE2n=-medianE2n/log(0.5);
    noisePower=meanE2n/EM_n;
    clear meanE2n;clear medianE2n;clear meanE2n;
    
 
    EstSumAn2=zero;
    EstSumAn2=EstSumAn2+ifftFilterArray.^2;
    
    EstSumAiAj=zero;
    EstSumAiAj=EstSumAiAj+ifftFilterArray;
    EstNoiseEnergy2=2*noisePower*sum(sum(EstSumAn2))+4*noisePower*sum(sum(EstSumAiAj));
    
    clear EstSumAn2;
    tau=sqrt(EstNoiseEnergy2/2);
    EstNoiseEnergy=tau*sqrt(pi/2);
    EstNoiseEnergySigma=sqrt((2-pi/2)*tau^2);    
    T=EstNoiseEnergy+k*EstNoiseEnergySigma;
    
    
    clear EstNoiseEnergy;clear EstNoiseEnergySigma;clear tau;
    clear EstNoiseEnergy2;clear EstSumAiAj;clear noisePower;
    
    T=T/1.7;
    Energy_ThisOrient=max(Energy_ThisOrient-T,zero);
    
    PC{1}=Energy_ThisOrient./(sumAn_ThisOrient+epsilon);  
  
    
 PC_temp =  PC{1};   
 for row = 1 :size(PC_temp,1)
    for line = 1 : size(PC_temp,2)
        if isnan(PC_temp(row,line))
           PC_temp(row,line) = double(0);
        end
    end
end
PC{1} =  PC_temp;

 for o=1:6   
    angl=(o-1)*pi/6;
    covx=PC{1}*cos(angl);
    covy=PC{1}*sin(angl);    
    %这是个累加操作
    covx2=covx2+covx.^2;
    covy2=covy2+covy.^2;
    covxy=covxy+covx.*covy;
 end
 
    denom=sqrt(covxy.^2+(covx2-covy2).^2)+epsilon;
    M=(covy2+covx2+denom)/2;
    m=(covy2+covx2-denom)/2;   
    %figure;imshow(M);   
    phaseCongruency=double(Energy_ThisOrient./(sumAn_ThisOrient+epsilon));   
    or=phaseCongruency;
    
for row = 1 :size(M,1)
    for line = 1 : size(M,2)
        if isnan(M(row,line))
            M(row,line) = double(0.00001);
            m(row,line) = double(0.00001);
            phaseCongruency(row,line) = double(0.00001);
            or(row,line) = double(0.00001);
        end
    end
end


%   figure;imshow(mapminmax(M,0,1));
%   figure;imshow(mapminmax(phaseCongruency,0,1));


end