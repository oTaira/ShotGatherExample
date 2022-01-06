 % create a shot gather with direct, head1, and head 2
 % part 1
clear; close all

i0 = 1:.1:89;
t1 = (2/50*log(cot(i0/2)));
for l=1:1:length(t1) %filter time to 500ms
    if (t1(l)<=.5 && t1(l) >= 0)
        time(l) = t1(l);
    end
end

x = (2*500/50*sinh(50*time/2));
z = (500/50*(cosh(50*time/2)-1));


N = 500; % number of time samples
numx = 48; % number of offsets
dx = 2; % offset spacing (m)
 
V1 = ones(1,881) .* 500; % velocity in 1st layer
V2 = 500+50*z; % velocity in 2nd layer
z12 = z; % depth of interface
% dt is the sampling interval (in seconds)
 
%Use Snell's law to find the critical angle
ic = asind(V1./V2);
%Calculate the critical distance
xc = 2*z12.*tand(ic);
 
 
dt = .002;
flag1 = 1; % 1=zero phase % 0 = minimum phase
i=[1:N];
t = (i-1)*dt;
final = zeros(N,1);
finaldelt = zeros(N,1);
trace = zeros(N, 1)';
%snr = 0.5; %extra credit
%noise = rand(N,1)' * 1/snr;
for k=1:numx
    x = (k-1)*dx;
    tdirect = x./V1;
    thead = 2.*z12.*cosd(ic)/V1 + (x-xc)/V2;
    treflect = (x.^2+4.*z12.^2).^(1/2)/V1;
    r = treflect .* V1;
    for j=1:N
        trace(j) = 0.;
        flag = abs(t(j)-tdirect) < dt/2;
        flagy = abs(t(j)-thead) < dt/2;
        flagr = abs(t(j)-treflect) < dt/2;
        if x<=xc % make head wave zero before critical distance
            flagh = 0;
        else
            flagh = 1;
        end
        if flag ==1
            trace(j) = 1.0;
        end
        if flagy == 1 && flagh==1
            trace(j) = 1.0;
        end
        if flagr == 1
            trace(j) = 1.0; %* 1/r; This is for question 6 
        end
        %trace(j) = trace(j) + noise(j); extra credit
        
    end    
if flag1 == 0
fmin = 10; %  wavelet low frequency (Hz)
fmax = 40;  %  wavelet high frequency (Hz)
[numer,denom]=butter(2,[fmin*dt*2,40*dt*2]);
result = filter(numer,denom,trace);
elseif flag1 == 1
fmin = 10; %  wavelet low frequency (Hz)
fmax = 75;  %  wavelet high frequency (Hz)
[numer,denom]=butter(2,[fmin*dt*2,fmax*dt*2]);
result = filtfilt(numer,denom,trace);
end
finaldelt = [finaldelt,trace'];
final = [final,result'];
end
final = final(:,2:49);
finaldelt = finaldelt(:,2:49);
wiggle(t, 0:dx:(numx-1)*dx, final)%, 'I'); This is when you dont want the 
hold on
%wiggle(t, 0:dx:(numx-1)*dx, finaldelt);
xlabel('Offset (m)')
ylabel('Time (s)') 
title(strcat('fmax =',  num2str(fmax), 'Hz'))

 %part 2
i0 = 1:.1:89;
x = real(2*500/50*sinh(50*time/2));
z = real(500/50*(cosh(50*time/2)-1));
xvalues = 10:.1:100;
xinterpt = real(interp1(x,xvalues));
tinterpt = real(asinh(xinterpt/1000)/25);
zinterpt = real(500/50*(cosh(50*tinterpt/2)-1));
figure(2)
plot(xinterpt, zinterpt);
 

%part 3

%increasing the velocity gradient decreases depth

%question 2
%part 1

for i=1:1:length(time)
    if z(i) >= 20
        V2(i) = 3000;
        
    end
    
end


fun = @(z) 1./((500 + 50.*z).*tan(500./(500 + 50.*z)));
treflection = integral(fun,0,20);
totalreflectiontime = treflection * 2;
criticalangle = asind(1500/3000);
totalrefractiontime = 2/50*log(cotd(criticalangle/2));

%part 2

i0 = 1:.1:89;
time = (2/50*log(cot(i0/2)));


x = (2*500/50*sinh(50*time/2));
z = (500/50*(cosh(50*time/2)-1));


N = 500; % number of time samples
numx = 48; % number of offsets
dx = 2; % offset spacing (m)
 
V1 = ones(1,881) .* 500; % velocity in 1st layer
V2; % velocity in 2nd layer
z12 = z; % depth of interface
% dt is the sampling interval (in seconds)
 
%Use Snell's law to find the critical angle
ic = asind(V1./V2);
%Calculate the critical distance
xc = 2*z12.*tand(ic);
 
 
dt = .002;
flag1 = 1; % 1=zero phase % 0 = minimum phase
i=[1:N];
t = (i-1)*dt;
final = zeros(N,1);
finaldelt = zeros(N,1);
trace = zeros(N, 1)';
%snr = 0.5; %extra credit
%noise = rand(N,1)' * 1/snr;
for k=1:numx
    x = (k-1)*dx;
    tdirect = x./V1;
    thead = 2.*z12.*cosd(ic)/V1 + (x-xc)/V2;
    treflect = (x.^2+4.*z12.^2).^(1/2)/V1;
    r = treflect .* V1;
    for j=1:N
        trace(j) = 0.;
        flag = abs(t(j)-tdirect) < dt/2;
        flagy = abs(t(j)-thead) < dt/2;
        flagr = abs(t(j)-treflect) < dt/2;
        if x<=xc % make head wave zero before critical distance
            flagh = 0;
        else
            flagh = 1;
        end
        if flag ==1
            trace(j) = 1.0;
        end
        if flagy == 1 && flagh==1
            trace(j) = 1.0;
        end
        if flagr == 1
            trace(j) = 1.0; %* 1/r; This is for question 6 
        end
        %trace(j) = trace(j) + noise(j); extra credit
        
    end    
if flag1 == 0
fmin = 10; %  wavelet low frequency (Hz)
fmax = 40;  %  wavelet high frequency (Hz)
[numer,denom]=butter(2,[fmin*dt*2,40*dt*2]);
result = filter(numer,denom,trace);
elseif flag1 == 1
fmin = 10; %  wavelet low frequency (Hz)
fmax = 75;  %  wavelet high frequency (Hz)
[numer,denom]=butter(2,[fmin*dt*2,fmax*dt*2]);
result = filtfilt(numer,denom,trace);
end
finaldelt = [finaldelt,trace'];
final = [final,result'];
end
final = final(:,2:49);
finaldelt = finaldelt(:,2:49);
figure(3)
wiggle(t, 0:dx:(numx-1)*dx, final)%, 'I'); This is when you dont want the 
hold on
%wiggle(t, 0:dx:(numx-1)*dx, finaldelt);
xlabel('Offset (m)')
ylabel('Time (s)') 
title(strcat('fmax =',  num2str(fmax), 'Hz'))
 
%part 3
i0 = 1:.1:89;
x = real(2*500/50*sinh(50*time/2));
xvalues = 10:.1:100;
xinterpt = real(interp1(x,xvalues));
tinterpt = real(asinh(xinterpt/1000)/25);
zinterpt = real(1500/50*(cosh(50*tinterpt/2)-1));
figure(5)
plot(xinterpt, zinterpt);

%question 3
%part 1
for l=1:1:length(t1) %filter time to 500ms
    if (t1(l)<=.5 && t1(l) >= 0)
        time(l) = t1(l);
    end
end
timetraveled = zeros(1,881);
for i=2:1:length(time)
timetraveled(i) = abs(time(i)-time(i-1));
end

vavg = sum(V2.*timetraveled) ./ sum(timetraveled);

%qustion 4
dx = 120;
dt = .1;
vs = dx/dt;
frequency = 1/(.08-.06);
