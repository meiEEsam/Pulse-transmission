close
clear
fc = 6.56e9;
Fs=5e10;
t=-50e-9:1/Fs:50e-9;
%%
%create pulse
pulse=zeros(1,length(t));
for i=2250:1:2750
    pulse(1,i)=1;
end
plot(t,pulse)
ylim([-1,2])
xlabel('time')
title('10 ns pulse')
%%      
%modulate the pulse
w0=10e9*2*pi;
mpulse=exp(1i*w0*t).*pulse;
%%
%take fft of modulated pulse
fftpulse=fft(mpulse);
L=length(pulse);
f = Fs*(0:floor((L))-1)/L;
plot(f,abs(fftpulse))
xlabel('frequency')
title('fft of the modulated pulse')
%%
%cancel frequancies below cutoff
for i=1:length(f)
    if (f(i)<fc)
        fftpulse(i)=0;

    end
end
figure
plot(f,abs(fftpulse))
%%
m=1;
for i=1:length(f)
    if (fftpulse(i)==0)
        m=m+1;
    else
        break
    end
end
f2=f(m:end);
fftpulse2=fftpulse(m:end);
%calculating propogation constants
beta=sqrt(1-(fc./f2).^2);
%shape of pulse in z=10
z=1000;
for i=1:length(fftpulse2)
    fftpulse2(i)=fftpulse2(i).*exp(-1i*z*beta(i));
end

%calculate ifft to find the result signal
fftpulse3=zeros(1,L);
fftpulse3(m:end)=fftpulse2;

figure
plot(f,abs(fftpulse3))

%%
%calculate ifft to find the result signal
disppulse=ifft(fftpulse3);
disppulse2=disppulse.*exp(-1i*w0*t);
plot(t,abs(real(disppulse2)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Gaussian pulse
%%
w0=10e9*2*pi;
t0=10e-9;
fc = 7.56e9;
Fs=5e10;
t=-150e-9:1/Fs:150e-9;
Gpulse=exp(-1*t.*t/(t0*t0));
plot(t,Gpulse)
ylim([-0.5,1.5])
title('Gaussian pulse')

%%

L=length(Gpulse);
f = Fs*(0:floor((L))-1)/L;
mGpulse=exp(1i*w0*t).*Gpulse;
plot(t,real(mGpulse))
ylim([-0.5,1.5])
title(' modulated Gaussian pulse')

%%
fftpulse=fft(mGpulse);
L=length(Gpulse);
f = Fs*(0:floor((L))-1)/L;
plot(f,abs(fftpulse))
xlabel('frequency')
title('fft of the modulated Gpulse')
close
%%
%cancel frequancies below cutoff
for i=1:length(f)
    if (f(i)<fc)
        fftpulse(i)=0;

    end
end
figure
plot(f,abs(fftpulse))
%%
m=1;
for i=1:length(f)
    if (fftpulse(i)==0)
        m=m+1;
    else
        break
    end
end
f2=f(m:end);
fftpulse2=fftpulse(m:end);
%calculating propogation constants
beta=sqrt(1-(fc./f2).^2);
%shape of pulse in z=10
z=15000;
for i=1:length(fftpulse2)
    fftpulse2(i)=fftpulse2(i).*exp(-1i*z*beta(i));
end

%calculate ifft to find the result signal
fftpulse3=zeros(1,L);
fftpulse3(m:end)=fftpulse2;
%%
%calculate ifft to find the result signal
disppulse=ifft(fftpulse3);
disppulse2=disppulse.*exp(-1i*w0*t);
%%
[a,b]=max(abs(real(disppulse2)));
[c,d]=max(abs(real(Gpulse)));
shiftpulse=zeros(1,L);
shiftpulse(d-b+1:end)=disppulse2(1:end-(d-b));
%shiftpulse(1:end-(b-d)+1)=disppulse2(b-d:end);
t0=0.65*10e-9;
Gpulse2=exp(-1*t.*t/(t0*t0));
plot(t,Gpulse2,'b')
hold on
plot(t,abs(real(shiftpulse)/1.1),'r')
title('z=15km');
legend('original pulse','output pulse')
ylim([-0.3,1.5])