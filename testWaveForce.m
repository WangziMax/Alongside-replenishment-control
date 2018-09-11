clc 
clear

i           = 0;
steps       = 720;
h           = 0.5;
w0          = 0.8;
x_w         = [0 0]';
Store_fw    = zeros(steps,1);
Store_fw(1) = 0;
U           = 10;
beta        = 45*pi/180;
fw2         = 1;
for i = 2:1:steps
	[fw1, x_w] = WaveForce(w0, h, x_w, beta, U);
    fw = fw1 + fw2;
	Store_fw(i)   = fw;
end

plot(Store_fw)