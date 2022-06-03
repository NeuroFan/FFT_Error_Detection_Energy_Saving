clear
clc

ERR_THRESHOLD = 2^-64; %anything larger than this is an error
N = 1024 ; %signal size


t = (0:N-1)*1/1e3;
y = chirp(t,0,1,250);

K = 1000;%1000
FaultRate = 0.0000000000000001;
errParseval = zeros(K,1); 
errActual = zeros(K,1); 
errHermitian = zeros(K,1); 
x = y';
f = fft(x); %really fft is not needed but just to make thing match

error_variance = zeros(K,1);
for i = 1:K
   disp(i)

  error_variance(i) =FaultRate;
  for j = 1:100

    %% error injection
    eR_indx = rand(N,1)<FaultRate;
    eI_indx = rand(N,1)<FaultRate;
    RandNoise =  random('Uniform',-1,1,[N,1]);
    eR = eR_indx .*RandNoise;
    eI = eI_indx .* RandNoise;
    f_injected = f + complex(eR,0) + complex(0,eI);
    %f_injected = f + complex(abs(real(f)).*eR,0) + complex(0,abs(imag(f)).*eI);  % the error is scaled to size of the bin so it is inpud dependant
    %%
    %eR = random('Uniform',0,FaultRate,[N,1]);
    %eI = random('Normal',0,FaultRate,[N,1]);
    %f_injected = f + complex(eR,0) + complex(0,eI);
    %% actual error measurement
    if abs(sum((f_injected - f).^2)./sum(f.^2))>ERR_THRESHOLD % if any error is injected at all 
         errActual(i)= errActual(i) + 1;
    end
    %% parseval's identity
    Ex = sum(x.^2);
    Ef = sum(abs(f_injected).^2)/N;
    ErrPercent = abs(Ex-Ef);
    if ErrPercent>ERR_THRESHOLD
        errParseval(i) = errParseval(i) + 1;
    end
  %% Hermitian error detection  
    R = real(f_injected);
    I = imag(f_injected);
     
    RH = real([0; f_injected((N:-1:2))]);
    IH = imag([0; f_injected((N:-1:2))]);
    
    ErReal = abs(R - RH);
    ErImag = abs(I + IH);
    
    if (max(ErReal(2:N))>ERR_THRESHOLD) || (max(ErImag(2:N))>ERR_THRESHOLD)   % if a Hermitian error is more than 1%
        errHermitian(i) = errHermitian(i) + 1;
    end
  end
    FaultRate = FaultRate+0.00001 ; %+0.00001
end
%% ploting
close all



subplot(2,1,1)
errActual_smoothed = movmean(errActual,100);
errParseval_smoothed = movmean(errParseval,100);
errHermitian_smoothed = movmean(errHermitian,100);


plot(error_variance,errActual_smoothed./j,'b')
hold on
plot(error_variance,errParseval_smoothed./j,'r')
hold on
plot(error_variance,errHermitian_smoothed./j,'g')


subplot(2,1,2)

ErrRate = errActual;
Parsevals_Silent_Error = errActual_smoothed - errParseval_smoothed;
Hermitians_Silent_Error =errActual_smoothed - errHermitian_smoothed;

Parsevals_Silent_Error_smoothed = movmean(Parsevals_Silent_Error./j,1);
Hermitians_Silent_Error_smoothed = movmean(Hermitians_Silent_Error./j,1);


plot(error_variance,Parsevals_Silent_Error_smoothed,'r');
hold on
plot(error_variance,Hermitians_Silent_Error_smoothed,'g');
