function [y_w, x_w] = WaveDisplacement(w0, h, x_w_last, beta, U)
% -----------------------------------------------------------------------------------------
% LINEAR WAVE MODEL based on 
% lambda: relative damping factor
% w0: the domninating wave frequency (rad/s)
% sigma: a constant describing the wave intensity
% sea state: 5 -> Hs = 3.5
% -----------------------------------------------------------------------------------------
Spectra = 'MPM';

%% Modified Pierson-Moskowitz spectra
if strcmp(Spectra, 'MPM')
	if w0 == 1.4
		lambda = 0.2606;
		sigma  = 0.8848;
	elseif w0 == 1.1
		lambda = 0.2588;
		sigma  = 0.9982;
	elseif w0 == 0.8;
		lambda = 0.2573;
		sigma  = 1.1705;
	end
	w0 = w0-(U*w0^2*cos(beta))/9.18;  % encounter frequency
end



Kw     = 2*lambda*w0*sigma;
                   
Aw 	   = [0			1
      	  -w0^2    	-2*lambda*w0];
Bw     = [0 Kw]';

x_w_dot = Aw*x_w_last + Bw*0.5*randn(1);

x_w = x_w_last + x_w_dot*h;

y_w = [0 1] * x_w;


%%
% Spectra = 'Tor';
% %% Torsethaugen spectra
% if strcmp(Spectra, 'Tor')
% 	if w0 == 1.4
% 		lambda = 0.1108;
% 		sigma  = 0.9932;
% 	elseif w0 == 1.1
% 		lambda = 0.1643;
% 		sigma  = 1.0808;
% 	elseif w0 == 0.8;
% 		lambda = 0.1850;
% 		sigma  = 1.3279;
% 	end
% end