%> @file allan_variance_testing.m
%> @brief test script to perform allan deviation analysis, mainly used to
%> validate first order Gauss-Markov IMU errors
%>
%> @author Eric Victorson
%> @date 2018/10/11
%> @version 1.0
% ======================================================================
%>
%> this testing script is used to verify the gauss markov noise created for
%> synthetic imu messages based on the __ error model for bias, bias
%> inrun, scale factor, and random walk
%>
%> when each of these processes are taken separately we should be able to
%> back out the markov noise amplitude (qc, or markov process standard deviation)
%> and correlation time (Tc, or time constant) according to the equations
%> c.16 - 18 and Figure C.6 of IEEE std 953-1997.
%>
%> educational sources:
%> https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=660628  IEEE std 953-1997
%> http://home.engineering.iastate.edu/~shermanp/AERE432/lectures/Rate%20Gyros/14-xvagne04.pdf
%> http://cache.freescale.com/files/sensors/doc/app_note/AN5087.pdf
%> https://openi.nlm.nih.gov/detailedresult.php?img=PMC3812568_sensors-13-09549f4&req=4
%> https://etd.ohiolink.edu/!etd.send_file?accession=osu1420709961&disposition=inline#page=144&zoom=100,0,661
%> https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=494457 IEEE std 647-1995 (bad qc and tc calculation)
%> note: when I say bad qc and tc calculation, I mean I believe IEEE to have published a mistake in this 
%> version of the standard.
% ======================================================================
%% Generate gauss markov processes to be analyzed

len = 1000000;
gyro_in = zeros(len, 3);
accel_in = zeros(len, 3);
dt = 1;

gmp1c = GMP1('IMU1');
[gyro_errors, accel_errors, s1_gyro_test, b1_gyro_test, b2_gyro_test, rw_gyro_test,...
                s1_accel_test, b1_accel_test, b2_accel_test, rw_accel_test] = gmp1c.run_imu(gyro_in, accel_in, dt);
            

accel = accel_errors(:,1);
bias = b2_accel_test(:,1);
bias_inrun = b1_accel_test(:,1);
scale = s1_accel_test(:,1);

%% Perform allan deviation analysis and create plots

fs = 1/dt;
pts = 200;

[T_bias,sigma_bias] = allan(bias,fs,pts);
[T_bias_inrun,sigma_bias_inrun] = allan(bias_inrun,fs,pts);
[T_scale,sigma_scale] = allan(scale,fs,pts);
[T_accel,sigma_accel] = allan(accel,fs,pts);

plot_allan_dev(T_bias, sigma_bias);
plot_allan_dev(T_bias_inrun, sigma_bias_inrun);
plot_allan_dev(T_scale, sigma_scale);
plot_allan_dev(T_accel, sigma_accel);


%% Calculate the qc tc terms
% first generate the theoretical allan variance for a gauss markov process
% with specified Tc and qc and plot along side the above created process

Tc = 300;
qc = 0.004899;
tau = linspace(0.1,10000,10000);
% generated gauss markov allan variance (exponentially correlated noise)
y_cn =  gen_avar_gauss_markov(qc, Tc, tau);
% numerical derivative of generated gauss-markov allan variance
d_y_cn = differentiate(tau, y_cn);

% plot the model and it's derivative, as well as the generated gauss markov
% process as well as its derivative
figure
loglog(tau, y_cn)
hold on
loglog(tau, d_y_cn)
loglog(T_bias_inrun, sigma_bias_inrun.^2)
loglog(T_bias_inrun, d_sigma_bias_inrun)
title('Gauss Markov Allan Variance')

grid on
legend('ycn', 'ycn numerical deriv', 'sigma bias inrun', 'd sigma bias inrun')

% estimate qc and Tc for the gauss markov process generated and for the
% allan variance definition

[qc_calc,tc_calc] = eval_correlated_noise_adev(T_bias_inrun, sigma_bias_inrun)
[qc_ycn,tc_ycn] = eval_correlated_noise_avar(tau, y_cn)
% note, for tau << Tc the equation for exponentially correlated markov noise reduces to the equation for rate random walk,
% so qc may be obtained by fitting a line to the +1/2 slope and finding
% where it intersects with tau = 3 (if only gauss markov noise is present)


%% provided allan deviation scripts
% fs = 1;
% [av, tau] = av_calc2(bias_inrun_tot, fs); % input rate in units deg/hr, fs in units Hz 
% av_fit_reduced_plotting(tau, av, 'sec', 'deg/hr');

%% plot sf, bias (b2), bias in run (b1) in normal data sheet units to compare
% to the output of the SECTR kalman filter estimates (likely not helpful)
dt = 0.01;
micro_g_to_g = 0.000001;
gravity_mean = 9.7976446561;
g_to_meters_per_sec2 = gravity_mean;
hours_to_seconds = 3600;
minutes_to_seconds = 60;
ppm_to_part = 0.000001;
deg_to_rad = pi / 180.0;
deg_per_hr_to_rad_per_sec = deg_to_rad / hours_to_seconds;
per_root_hour_to_per_root_second = 1.0/60.0;

x = linspace(0,length(s1_gyro)*0.01, length(s1_gyro));
s1_gyro_p = s1_gyro .*(1/ppm_to_part);
s1_accel_p = s1_accel .*(1/ppm_to_part);
b1_gyro_p = b1_gyro .*(1/dt) .*(1/deg_per_hr_to_rad_per_sec);
b2_gyro_p = b2_gyro .*(1/dt) .*(1/deg_per_hr_to_rad_per_sec);
b1_accel_p = b1_accel .*(1/dt) .*((1/g_to_meters_per_sec2)*(1/micro_g_to_g));
b2_accel_p = b2_accel .*(1/dt) .*((1/g_to_meters_per_sec2)*(1/micro_g_to_g));


figure
% gyroscope
subplot(2,2,1)
plot(x, b2_gyro_p)
title('gyro bias');
legend('x', 'y', 'z')
xlabel('Time (s)')
ylabel('Degrees / hr')
grid on

subplot(2,2,2)
plot(x, b1_gyro_p)
title('gyro inrun bias');
legend('x', 'y', 'z')
xlabel('Time (s)')
ylabel('Degrees / hr')
grid on

subplot(2,2,3)
plot(x, s1_gyro_p)
title('gyro scale factor');
legend('x', 'y', 'z')
xlabel('Time (s)')
ylabel('PPM')
grid on

% accelerometer
figure

subplot(2,2,1)
plot(x, b2_accel_p)
title('accel bias');
legend('x', 'y', 'z')
xlabel('Time (s)')
ylabel('Micro G')
grid on

subplot(2,2,2)
plot(x, b1_accel_p)
title('accel inrun bias');
legend('x', 'y', 'z')
xlabel('Time (s)')
ylabel('Micro G')
grid on

subplot(2,2,3)
plot(x, s1_accel_p)
title('accel scale factor');
legend('x', 'y', 'z')
xlabel('Time (s)')
ylabel('PPM')
grid on

%% accel bias inrun autocorrelation
% the autocorrelation plot can easily be used to verify the correlation time,
% Tc.  Tc will occur when a long enough enough sample is run (long enough
% for the correlated terms to be evident), at the value that causes the
% autocorrelation to be 1/e its initial value.  The variance that is used
% to define the autocorrelation function is DIFFERENT from the variance
% that defines the driven noise process.  The autocorrelation variance is
% the variance for the entire time series.

% autocorr(b1_accel_p(:,1))
t = (0:length(b1_accel(:,1))-1)*dt;
T = 300.0;
bias_inrun_stddev = 300 * micro_g_to_g * g_to_meters_per_sec2;
% [acor, lag] = xcorr(b1_accel(:,1),b1_accel(:,1), 'coeff');
[acor, lag] = xcorr(b1_accel(:,1),b1_accel(:,1));

tau = linspace(0,450,450);
acor_expected = (bias_inrun_stddev^2) * exp(-1*(tau./T));
% acor = acor .* acor_expected(1); % what if I scale it?


figure
plot(lag*dt, acor)
hold on
plot(tau, acor_expected)
grid on

figure 
plot(b1_accel(:,1))


%% helper functions

% ======================================================================
%> @brief plot_allan plot the allan deviation
%> 
%> Plots the allan deviation as a function of sample cluster time.
%> performs line fitting and calculation of random walk, rate random walk,
%> bias instability.
%>
%> @param T sample cluster time vector, output from allan()
%> @param sigma allan deviation vector, output from allan()
% ======================================================================
function plot_allan_dev(T,sigma)

    c = @cmu.colors; % shortcut function handle
    deriv = differentiate_log(T, sigma);
    
    % find where slope is nearest to -0.5 (angle / velocity random walk)
    % (intersection with line at tau = 1)
    [v1,idx1] = find_nearest(deriv, -0.5);
    
    % find where slope is nearest to +0.5 (rate random walk) (intersection with
    % line at tau = 3)
    [val2,idx2] = find_nearest(deriv, 0.5);
    
    % RATE RANDOM WALK (RED NOISE (white noise accumulation) f^-2)
    % asymptotic properties of allan variance
    % sigma^2 = ((K^2)/3)*tau
    % 
    % Long term changes to bias offset will be randomly distributed and may be permanent in nature. Even
    % though the drift of an individual sensor can not be predicted, the time scale over which the changes
    % occur can be defined by the RRW and introduces the opportunity to plan for recalibration in critical
    % applications that require extended life
    % This equation shows that RRW is represented by +0.5 slope on a log-log plot of ?(? ) (fig. 3). RRW
    % constant K can be read off from plot at ? = 3.
    A1_rrw = log10(sigma(idx2)/(T(idx2)^0.5));
    
    % rate random walk is caused by the bias changing slowly with time
    RRW1_line = 10^(A1_rrw)*T.^0.5;

    RRW1 = 10^(A1_rrw)*3.^0.5;
    
    % VELOCITY RANDOM WALK / ANGLE RANDOM WALK (WHITE NOISE (johnson nyquites thermal noise) f^0)
    % asymptotic properties of allan variance
    % sigma^2 = (N^2) / tau
    %
    % ARW is high frequency noise and it can be observed as the short-term variation in the output. After
    % performing an integration, it causes random error in angle with distribution, which is proportional to
    % the square root of the elapsed time
    % 
    % The equation (9) shows that the slope of ?(? ) log-log plot is ?0.5 and the coefficient N can be
    % obtained from the plot at ? = 1
    % 
    % White noise (also referred to as angle random walk [25]),
    % originating from, e.g., thermal noise [27], is a simple stochastic
    % process that consists of independent and identically distributed
    % samples. Its variance is inversely proportional to the averaging
    % time, thus the effect of white noise can be s

    % why doesn't this work when I look at comparing derivatives?
    %idx1 = 1;
    %idx3 = 1;

    A1 = log10(sigma(idx1)/(T(idx1)^-0.5));

    VRW_line = 10^(A1)*T.^-0.5;

    VRW = 10^(A1)*1.^-0.5;
    
    % BIAS INSTABILITY (PINK NOISE (ELECTRONICS FLICKER) f^-1)
    % asymptotic properties of allan variance
    % ~ 2*(B^2)*ln(2)/pi
    %
    % The bias instability has an impact on the long-term stability. It is slow fluctuation of output so it
    % appears in low frequencies as 1/f noise. Bias instability determines the best stability that could be
    % achieved with fully modeled sensor and active bias estimation.  Van be
    % used to measure the power of 1/f noise

    [accel_bias_instability,idx_accel] = min(sigma);
    
%     accel_VRW_text = text(1, VRW, num2str(VRW));
% 
%     accel_bias_instability_text = text(T(idx_accel), accel_bias_instability, num2str(accel_bias_instability));

    figure
    loglog(T,sigma, 'color',c('deep carrot orange') )
    hold on
    loglog(T, VRW_line, '--', 'color', c('medium electric blue'))
    loglog(1,VRW,'s', 'MarkerEdgeColor',c('medium electric blue'),'MarkerFaceColor',c('medium electric blue'))
    loglog(T(idx_accel), accel_bias_instability,'s', 'MarkerEdgeColor',c('deep carrot orange'),'MarkerFaceColor',c('deep carrot orange'));
    loglog(T, RRW1_line,'--', 'color', c('jungle green'));
    loglog(3,RRW1,'s', 'MarkerEdgeColor',c('jungle green'),'MarkerFaceColor',c('jungle green'))
    grid on
    text(T(idx_accel), accel_bias_instability*1.4, num2str(accel_bias_instability));
    text(1.1, VRW*1.3, num2str(VRW));
    text(3, RRW1*0.7, num2str(RRW1));
    title('Accel Allan Deviation vs Sample Period');
    xlabel('Sample Period (s)');
    ylabel('Allan Deviation');
    legend('Allan Deviation', 'RW fit line', 'Random Walk (m/s/s/100hz)', 'Bias Instability (m/s/s)', 'RRW fit line', 'RRW (m/s/s*sqrt(100hz))')
    hold off
end

% ======================================================================
%> @brief allan calculate allan deviation
%> 
%> Generates the allan deviation of the input time series, adapted from:
%> http://cache.freescale.com/files/sensors/doc/app_note/AN5087.pdf
%>
%> @param omega input time series vector (either gyro rate output, or
%accelerometer rate output)
%> @param fs sampling frequency
%> @param pts number of points to use for log spaced averaging factor
%> vector.  pts can be chosen arbitrarily such that pts < (#samples - 1)/2
%>
%> @retval T allan deviation averaging factor
%> @retval sigma allan deviation
% ======================================================================
function [T,sigma] = allan(omega,fs,pts)
    [N,M] = size(omega);                                % figure out how big the output data set is
    n = 2.^(0:floor(log2(N/2)))';                       % determine largest bin size
    maxN = n(end);
    endLogInc = log10(maxN);
    
    m = unique(ceil(logspace(-3,endLogInc,pts)))';      % create log spaced vector average factor
    t0 = 1/fs;                                          % t0 = sample interval
    T = m*t0;                                           % T = length of time for each cluster

    theta = cumsum(omega)/fs;                           % cumulative sum
    sigma2 = zeros(length(T),M);                        % array of dimensions (cluster periods) X (#variables)
    
    for i=1:length(m)                                   % loop over the various cluster sizes
        for k=1:N-2*m(i)                                % implements the summation in the AV equation
            sigma2(i,:) = sigma2(i,:) + (theta(k+2*m(i),:) - 2*theta(k+m(i),:) + theta(k,:)).^2;
        end
    end

    sigma2 = sigma2./repmat((2*T.^2.*(N-2*m)),1,M);     % allan variance
    sigma = sqrt(sigma2);                               % allan deviation
end

% ======================================================================
%> @brief differentiate method to differentiate time series
%>
%> @param T time series to differentiate with respect to
%> @param sigma time series to be differentiated
%>
%> @retval d returns differentiated time series
% ======================================================================
function [d] = differentiate(T, sigma)
    for i = 1:(length(T)-1)
       d(i,1) = (sigma(i+1) - sigma(i)) / (T(i+1) - T(i));
    end
    d(end+1) = d(end);
end

% ======================================================================
%> @brief differentiate_log method to differentiate log scaled time series
%>
%> @param T time series to differentiate with respect to
%> @param sigma time series to be differentiated
%>
%> @retval d returns differentiated time series
% ======================================================================
function [d] = differentiate_log(T, sigma)
    T = log10(T);
    sigma = log10(sigma);
    for i = 1:(length(T)-1)
       d(i,1) = (sigma(i+1) - sigma(i)) / (T(i+1) - T(i));
    end
    d(end+1) = d(end);
end

function [val,idx] = find_nearest(in, value)
    [val, idx] = min(abs(in-value));
end

function out = add_white_noise(in, sigma)
    out = zeros(length(in),0);
    for i = 1 : length(in)
       out(i) = in(i) + normrnd(0,sigma); 
    end
end


% ======================================================================
%> @brief eval_correlated_noise_adev method to estimate qc and Tc for first
%> order gauss markov process from the allan deviation.
%>
%> NOTE: The IEEE standard for ring laser gyros (IEEE std 647-1995 / 2006) has an error in the
%> calculation of qc and Tc, instead use the IEEE standard for fiber optic
%> gyros (IEEE std 952-1997)
%> if only exponentially correlated noise is present, the maxima of the
%> allan deviation will occur when tau / Tc = 1.8926 (emperically proven ALWAYS).  This can
%> then be plugged into the overall gauss markov allan variance equation to
%> solve for sigma_max when tau/Tc = 1.89.  The equation then reduces to
%> sigma_max ~ 0.437 * qc * sqrt(Tc)
%> The relationship shown in figure C.6 of IEEE std 647-1995 / 2006 doesn't
%> make sense dimensionally
%>
%> @param T allan deviation sampling factor
%> @param sigma allan deviation
%>
%> @retval qc correlated noise estimated standard deviation
%> @retval Tc correlated noise estimated time constant
% ======================================================================
function [qc, Tc] = eval_correlated_noise_adev(T, sigma)
    sigma2 = sigma.^2;
    [maxs, maxs_idx] = max(sigma2);
    Tc = T(maxs_idx) / 1.8926;
    qc = sqrt(sigma2(maxs_idx)) / (0.4368 * sqrt(Tc));
end

% ======================================================================
%> @brief eval_correlated_noise_avar method to estimate qc and Tc for first
%> order gauss markov process from the allan variance
%> @param T allan deviation sampling factor
%> @param sigma allan deviation
%>
%> @retval qc correlated noise estimated standard deviation
%> @retval Tc correlated noise estimated time constant
% ======================================================================
function [qc, Tc] = eval_correlated_noise_avar(T, sigma2)
    [maxs, maxs_idx] = max(sigma2);
    Tc = T(maxs_idx) / 1.8926;
    qc = sqrt(sigma2(maxs_idx)) / (0.4368 * sqrt(Tc));
end

% ======================================================================
%> @brief gen_avar_gauss_markov generate expected allan variance for gauss
%> markov process
%> 
%> @param qc gauss-markov driven noise amplitude (standard deviation)
%> @param Tc gauss-markov time constant
%> @param tau allan deviation sampling factor vector
%>
%> @retval y_cn expected allan variance time series
% ======================================================================
function y_cn = gen_avar_gauss_markov(qc, Tc, tau)
    y_cn =  (((qc.*Tc).^2)./tau).*(1.0 - ((Tc./2.0./tau).*(3.0 - 4.0*exp(-tau./Tc)+exp(-2.0.*tau./Tc))));
end
