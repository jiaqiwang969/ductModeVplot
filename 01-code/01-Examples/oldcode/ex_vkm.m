      clc
      clear
      fs = 4000;
      T = 5;
      dt = 1/fs;
      t = (0:dt:(T-dt))';
      N = numel(t);

      % Instationary component
      A1 = [linspace(0.5,1,floor(N/2)) linspace(1,0.5,ceil(N/2))]';
      f1 = 0.2*fs - 0.1*fs*cos(pi*t/T);
      phi1 = 2*pi*cumsum(f1)*dt;
      y1 = A1.*cos(phi1);

      % Stationary component
      A2 = ones(N,1);
      f2 = 0.2*fs*ones(N,1);
      phi2 = 2*pi*cumsum(f2)*dt;
      y2 = A2.*sin(phi2);

      % White noise
      e = 2*rand(size(y1));

      % Mixed signal
      y = y1 + y2 + e;

      % Perform VKF on periodic components
      p = 2;
      bw = 1;
      [a,c] = vkf(y,fs,[f1 f2],p,bw);
      x = real(a.*c);

      % Reveal white noise
      w = y-sum(x,2);

      % Plot
      figure('color','white')
      subplot(2,1,1)
      spectrogram(y,round(fs/16),[],[],fs)
      subplot(2,1,2), hold on
      plot(t,[A1 A2],'--')
      plot(t,abs(a))
      xlabel('Time (s)')
      ylabel('Amplitude')

      % Playback
      soundsc(y,fs) % Original signal
      soundsc(x,fs) % Periodic components (2 channels)
      soundsc(w,fs) % Remaining white noise