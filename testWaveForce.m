clc 
clear

i           = 0;
steps       = 720;
h           = 0.1;					% time step (s)
w0          = 1.1;					% wave frequency (rad/s)
x_w         = [0 0]';
Store_dw    = zeros(steps,1);
Store_dw(1) = 0;
Store_fw    = zeros(steps,1);
Store_fw(1) = 0;
U           = 10;					% speed
beta        = 45*pi/180;			% encounter angle

k1          = 0.6;
k2          = 0.8;
k3          = 0.3;

fw2         = [-1, -1, 0]';		% 2nd-order wave drift forces (N,N*m)

for i = 2:1:steps
	% the oscillatoric motion from 1st-order wave disturbances
	[disw1, x_w] = WaveDisplacement(w0, h, x_w, beta, U);
	Disw1        = [k1, k2, k3]'*disw1;							% oscillatoric motion
	Store_dw(i)  = disw1;

	% 2nd-order wave drift
	fw2_dot     = 0.1*randn(1)*ones(3,1);	
	fw2         = fw2 + h*fw2_dot;
	Store_fw(i) = fw2(1);
end
plot(Store_dw)
% plot(Store_fw)