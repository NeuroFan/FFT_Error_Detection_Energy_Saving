clear
clc

ERR_THRESHOLD = 1/10^9; %anything larger than this is an error
N = 1024 ; %signal size


t = (0:N-1)*1/1e3;
y = chirp(t,0,1,250);

K = 10000;%1000
FaultRate = 0.0000;
errParseval = zeros(K,1); 
errActual = zeros(K,1); 
errHermitian = zeros(K,1); 
x = y';
f = fft(x); %really fft is not needed but just to make thing match
 
for i = 1:K
   disp(i)
  FaultRate = FaultRate+0.000001 ; %+0.00001
  for j = 1:1000

    %% error injection
    eR_indx = rand(N,1)<FaultRate;
    eI_indx = rand(N,1)<FaultRate;
    RandNoise =  random('Normal',0,0.5,[N,1]);
    eR = eR_indx .*RandNoise;
    eI = i*eI_indx .* RandNoise;
    f_injected = f + abs(real(f)).*eR + abs(imag(f)).*eI;  % the error is scaled to size of the bin so it is inpud dependant

    %% actual error measurement
    if sum(abs(eR) + abs(eI))>0 % if any error is injected at all 
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
    RH = real([0; f_injected((N:-1:2))]);
    
    I = real(f_injected);
    IH = real([0; f_injected((N:-1:2))]);
    
    ErReal = abs(R - RH);
    ErImag = abs(I - IH);
    
    if (max(ErReal(2:N))>ERR_THRESHOLD) || (max(ErImag(2:N))>ERR_THRESHOLD)   % if a Hermitian error is more than 1%
        errHermitian(i) = errHermitian(i) + 1;
    end
  end
end
%% ploting
close all



subplot(2,1,1)
errActual_smoothed = movmean(errActual,100);
errParseval_smoothed = movmean(errParseval,100);
errHermitian_smoothed = movmean(errHermitian,100);


plot(errActual_smoothed./j,'b')
hold on
plot(errParseval_smoothed./j,'r')
hold on
plot(errHermitian_smoothed./j,'g')


subplot(2,1,2)

ErrRate = errActual;
Parsevals_Silent_Error = errActual_smoothed - errParseval_smoothed;
Hermitians_Silent_Error =errActual_smoothed - errHermitian_smoothed;

Parsevals_Silent_Error_smoothed = movmean(Parsevals_Silent_Error./j,1);
Hermitians_Silent_Error_smoothed = movmean(Hermitians_Silent_Error./j,1);


plot(Parsevals_Silent_Error_smoothed,'r');
hold on
plot(Hermitians_Silent_Error_smoothed,'g');
