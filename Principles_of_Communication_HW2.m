%%
% Name:         Bithiah Ngan
% Course:       ELC 4350 
% Assignment:   HW 2 
% Date:         3/20/2021
%%
clear all
clc

%TRANSMITTER
% encode text string as T-spaced PAM (+/-1, +/-3) sequence
str='01234 I wish I were an Oscar Meyer wiener 56789';
m=letters2pam(str); N=length(m);    % 4-level signal of length N
% zero pad T-spaced symbol sequence to create upsampled T/M-spaced
% sequence of scaled T-spaced pulses (with T = 1 time unit)
M=100; mup=zeros(1,N*M); mup(1:M:end)=m; % oversampling factor
% Hamming pulse filter with T/M-spaced impulse response
p=hamming(M);                       % blip pulse of width M
x=filter(p,1,mup);                  % convolve pulse shape with data
figure(1), plotspec(x,1/M)          % baseband signal spectrum
% am modulation
t=1/M:1/M:length(x)/M;              % T/M-spaced time vector
fc= 0.7;                              % carrier frequency
c=cos(2*pi*fc*t);                   % carrier
r=c.*x;                             % modulate message with carrier

%RECEIVER
% am demodulation of received signal sequence r
c2=cos(2*pi*fc*t);                   % synchronized cosine for mixing
x2=r.*c2;                            % demod received signal
fl=50; fbe=[0 0.1 0.2 1]; damps=[1 1 0 0 ];  % design of LPF parameters
b=firpm(fl,fbe,damps);               % create LPF impulse response
x3=2*filter(b,1,x2);                 % LPF and scale downconverted signal
% extract upsampled pulses using correlation implemented as a convolving filter
y=filter(fliplr(p)/(pow(p)*M),1,x3); % filter rec'd sig with pulse; normalize
% set delay to first symbol-sample and increment by M
z=y(0.5*fl+M:M:end);                 % downsample to symbol rate
figure(2), plot([1:length(z)],z,'.') % soft decisions
title('Soft Decision')
% decision device and symbol matching performance assessment
mprime=quantalph(z,[-3,-1,1,3])';    % quantize to +/-1 and +/-3 alphabet
cvar=(mprime-z)*(mprime-z)'/length(mprime), % cluster variance
lmp=length(mprime);
pererr=100*sum(abs(sign(mprime-m(1:lmp))))/lmp, % symb err
% decode decision device output to text string
reconstructed_message=pam2letters(mprime)      % reconstruct message

%%Functions
% f = pam2letters(seq)
% reconstruct string of +/-1 +/-3 into letters
function f = pam2letters(seq)

S = length(seq);
off = mod(S,4);

if off ~= 0
  sprintf('dropping last %i PAM symbols',off)
  seq = seq(1:S-off);
end

N = length(seq)/4;
f=[];
for k = 0:N-1
  f(k+1) = base2dec(char((seq(4*k+1:4*k+4)+99)/2),4);
end

f = char(f);
end

% f = letters2pam(str)
% encode a string of ASCII text into +/-1, +/-3

function f = letters2pam(str);           % call as Matlab function
N=length(str);                           % length of string
f=zeros(1,4*N);                          % store 4-PAM coding here
for k=0:N-1                              % change to "base 4"
  f(4*k+1:4*k+4)=2*(dec2base(double(str(k+1)),4,4))-99;
end
end

% plotspec(x,Ts) plots the spectrum of the signal x
% Ts = time (in seconds) between adjacent samples in x
function plotspec(x,Ts)
N=length(x);                               % length of the signal x
t=Ts*(1:N);                                % define a time vector
ssf=(ceil(-N/2):ceil(N/2)-1)/(Ts*N);       % frequency vector
fx=fft(x(1:N));                            % do DFT/FFT
fxs=fftshift(fx);                          % shift it for plotting
subplot(2,1,1), plot(t,x)                  % plot the waveform
title('The behavior of the ideal system at fc = 50')
xlabel('seconds'); ylabel('amplitude')     % label the axes
subplot(2,1,2), plot(ssf,abs(fxs))         % plot magnitude spectrum
xlabel('frequency'); ylabel('magnitude')   % label the axes

end
% y=pow(x) calculates the power in the input sequence x
function y=pow(x)
y=x(:)'*x(:)/length(x);
end

function y=quantalph(x,alphabet)
alphabet=alphabet(:);
x=x(:);
alpha=alphabet(:,ones(size(x)))';
dist=(x(:,ones(size(alphabet)))-alpha).^2;
[v,i]=min(dist,[],2);
y=alphabet(i);
end