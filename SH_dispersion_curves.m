%%
%%% SH waves dispersion curves based on formulas from Ultrasonic Waves in
%%% Solid Media by Joseph Rose
%close all
ct=3.23e3; % Shear velocity of Carbon Steel
%d=7.1628*1e-3; % Thickness of the plate
prompt = 'Thickness of plate[mm]: ';
d = input(prompt)*1e-3;

prompt = 'Rows of magnets: ';
n = input(prompt);

prompt = 'Magnet width [in]: ';
magw = (25.4*input(prompt)) + 0.125;
coil_length = (magw*n/25.4) %in inches

fmax=1000e3; % Maximum frequency for the dispersion curves
fd=linspace(0.00001,10e6*d,10000);
Lfd=length(fd);

nmax=10;
fd0=fd;
phC0=ct*ones(1,Lfd);
gC0=ct*ones(1,Lfd);
fdf=[fd0];
phC=[phC0];
gC=[gC0];
for n=1:nmax
    cof=n*ct/2;
    if cof<fd(Lfd)
        i=1;
        while fd(i)<(cof)
            i=i+1;
        end
        eval(strcat('fd',num2str(n),'=fd(i:Lfd);'));
        eval(strcat('phC',num2str(n),'=2*ct*fd',num2str(n),'./sqrt(4*fd',num2str(n),'.*fd',num2str(n),'-n*n*ct*ct);'));
        eval(strcat('gC',num2str(n),'=ct*sqrt(1-((n/2)^2)./((fd',num2str(n),'/ct).^2));'));
        eval(strcat('fdf=[fdf fd',num2str(n),'];'));
        eval(strcat('phC=[phC phC',num2str(n),'];'));
        eval(strcat('gC=[gC gC',num2str(n),'];'));
    else
        nmax=n;
    end
end

figure(1)
set(1,'position',[80 500 700 430],'Color',[1,1,1]);
plot(fdf/d*1e-3,phC,'k*');
xlabel('Frequency (kHz)');ylabel('Phase velocity (m/s)');
axis([fdf(1)/d*1e-3 fmax*1e-3 2500 12000])

%{
figure(2)
set(2,'position',[880 500 700 430],'Color',[1,1,1]);
plot(fdf/d*1e-3,gC,'k*');
xlabel('Frequency (kHz)');ylabel('Group velocity (m/s)');
axis([fdf(1)/d*1e-3 fmax*1e-3 0 4000])
%}
%% Plot Spacial Bandwidth
%n = 4;
%magw = 3.175+0.125;
T = 2*magw;
dx = magw/1000;
w = 2*pi/T;
x = 0:dx:n*magw;
F = sin(w*x);
Fsq = F.^2;

F1 = [F zeros(1,500000)];
L = length(F1);
NFFT = 2^nextpow2(L);
Y = fft(F1,NFFT);
Fs = 1/dx;
f = Fs/2*linspace(0,1,NFFT/2+1);

figure(1)
hold on
plot(fdf/d*1e-3,0.001*T*fdf/d,'k:');

%{
figure(3)
plot(f,abs(Y(1:NFFT/2+1))./(max(abs(Y(1:NFFT/2+1))))) 
hold on
plot(f,(abs(Y(1:NFFT/2+1))./(max(abs(Y(1:NFFT/2+1))))).^2,'r:') 
xlim([0 0.4])
xlabel('Spacial Frequency (cyc/mm)')
ylabel('Amplitude |.|')
legend('Response','Squared Response')
%}
%Find -6dB Bandwidth
amp = (abs(Y(1:NFFT/2+1))./(max(abs(Y(1:NFFT/2+1))))).^2;
[p,l] = max(amp);
s = length(amp);

bw6db = 10^(-6/10);
l1 = max(find(amp(1:l)<= bw6db));
l2 = min(find(amp(1:l)>= bw6db));
u1 = min(find(amp(l:s)<= bw6db))+l-1;
u2 = max(find(amp(l:s)>= bw6db))+l-1;
bw1=interp1(amp([l1 l2]),f([l1,l2]),bw6db);
bw2=interp1(amp([u1 u2]),f([u1,u2]),bw6db);

lambda6db1 = 1/bw1
lambda6db2 = 1/bw2

figure(1)
%plot(fdf/d*1e-3,0.001*lambda6db1*fdf/d,'r:');
%plot(fdf/d*1e-3,0.001*lambda6db2*fdf/d,'r:');

%% Time Dependant Frequency

prompt = 'Optimal frequency? ';
f0 = input(prompt)
lambda = T/1000;
E = n/2;
cg = 5000;
N = 5;
w0 = 2*pi*f0;
dt = (N/f0)/1000;
t = 0:dt:N/f0;
pulse = 34*sin(w0*t);

plot(t,pulse)

Fs = 1/dt;
P1 = [pulse zeros(1,500000)];
L = length(P1);
NFFT = 2^nextpow2(L);
X = fft(P1,NFFT);
f = Fs/2*linspace(0,1,NFFT/2+1);

%figure(4)
%plot(f,abs(X(1:NFFT/2+1))./(max(abs(X(1:NFFT/2+1))))) 

%xlim([200000 1600000])
%xlabel('Frequency (Hz)')
%ylabel('Amplitude |.|')


%Find -6dB Bandwidth
amp = (abs(X(1:NFFT/2+1))./(max(abs(X(1:NFFT/2+1))))).^2;
[p,l] = max(amp);
bw6db = 10^(-6/10);
s = length(amp);
l1 = max(find(amp(1:l)<= bw6db));
l2 = min(find(amp(1:l)>= bw6db));
u1 = min(find(amp(l:s)<= bw6db))+l-1;
u2 = max(find(amp(l:s)>= bw6db))+l-1;
bw1=interp1(amp([l1 l2]),f([l1,l2]),bw6db);
bw2=interp1(amp([u1 u2]),f([u1,u2]),bw6db);

lambda6db1 = bw1
lambda6db2 = bw2


%Find Main Lobe Width - Estimate Noise Floor To be 10%
nf = 0.1;
l1 = max(find(amp(1:l)<= nf));
l2 = min(find(amp(1:l)>= nf));
u1 = min(find(amp(l:s)<= nf))+l-1;
u2 = max(find(amp(l:s)>= nf))+l-1;

bw1=interp1(amp([l1 l2]),f([l1,l2]),nf);
bw2=interp1(amp([u1 u2]),f([u1,u2]),nf);

