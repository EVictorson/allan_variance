%> @file GMP1.m
%> @brief Class for modeling discrete time gauss markov errors for IMUs.
%>
%> @author Eric Victorson
%> @date 2018/10/11
%> @version 1.0
% ======================================================================
%> @brief GMP1 Class used to generate a first order gauss markov error process
%> to model IMU errors
%> Additional Details:
%> 
%> From the navigation class part 2:
%> Nav_class part 2 section 5-9
%> As previously stated, the gyro bias error is typically modeled as a "slow" 
%> Gauss-Markov process (1000 hour correlation time) which can be considered 
%> a random constant for a single flight.  The typical gyro bias error budget
%> is on the order of 0.01 deg/hr.
%> 
%> Nav class part 2 section 5-11
%> As stated before, the gyro scale factor error is also typically modeled as
%> a Gauss-Markov process with a long correlation time (1000 hours).  Again, 
%> for a single flight this can be considered as a random constant. 
%> The typical gyro scale factor error budget is on the order of 5 
%> parts-per-million RMS
%>
%> Again, as with the gyro, the accel bias error is typically modeled as a 
%> "slow" Gauss-Markov process which can be considered a random constant for 
%> a single flight.  The typical accel bias error budget is on the order 
%> of 84 µg.
%> 
%> The accel scale factor error is also typically modeled as a Gauss-Markov 
%> process with a long correlation time (100 hours).  Again, for a single 
%> flight this can be considered as a random constant. 
%> The typical accel scale factor error budget is on the order of 300 
%> parts-per-million RMS.
%>
%> The following terms are neglected, but documented for future use if
%> necessary: scale factor nonlinearity, nonorthogonality, misalignment
%>
%> Validation of this class has been performed via allan deviation
%> analysis, wherein each individual process (bias, bias inrun, scale factor)
%> was analyzed and the time constant (Tc) and noise magnitude (qc) were
%> verified.  For performing this type of analysis again see the
%> allan_variance_testing.m file.
%>
%> At it's core, an exponentially correlated, noise driven stochastic
%> process (gauss markov process) may be modeled with the following
%> differential equation: dx/dt = (-1/tau) * x + N(mu, sigma)
%> which is just an exponentially decaying autoregressive term with
%> driven white noise.
%> The discrete time version of this is: 
%> x(n+1) = x(n) * exp((-1/tau) * dt) + N(mu,sigma) * dt
%> 
%> The first order gauss-markov process may also be defined in terms of
%> its autocorrelation function: Rxx = (sigma^2)*exp(-beta*tau), where sigma
%> is the ENTIRE process variation (not the driven noise variation), beta is
%> the inverse of the correlation time T = 1/beta (time constant), and tau is
%> the autocorrelation time lag.  For sufficiently long time series (much 
%> longer than the time constant), the autocorrelation plot of the first order 
%> Gauss-Markov process can verify the correlation time, which is the time 
%> lag at which Rxx = (sigma^2) / e, and the entire time series variance 
%> will appear as the peak at time lag = 0.
%> Note that the driven noise standard deviation cannot be validated via
%> autocorrelation. Also note that matlab has two different autocorrelation
%> functions that will produce slightly different results if your time scales
%> are short (autocorr and xcorr), as well as an autocovariance function, xcov.
%> Autocorr and xcov are normalized by the process mean such that a constant
%> process will have a steady autocorrelation (autocorrelation is just
%> autocovariance scaled by the inverse of the process variance). xcorr,
%> however, does not do this, so a constant process will have a decaying
%> autocorrelation as the time lag increases, which results in a shortening
%> of the summation by 1 element each iteration.
%> 
%> For process lengths on the order of magnitude of the time constant or
%> less the process will be dominated by the integrated white noise term, resulting
%> in a random walk.  For time scales much longer than the time constant the
%> exponentially correlated term will begin to be apparent.  This can be see
%> on an allan deviation plot, where for sampling intervals much shorter than
%> the time constant the gauss-markov allan variance reduces to that of a
%> singly integrated white noise process (rate random walk), whose slope 
%> is +1/2, and the noise magnitude (standard deviation) may be picked off 
%> by finding the intersection of the +1/2 slope line and sampling_interal = 3.  
%> This can also be verified by differentiating this time series to obtain 
%> gaussian noise, whose allan deviation has a slope of -1/2, and the noise
%> magnitude is interpreted as the intersection of thise line with
%> sampling_interval = 1.
%>
%> For process time scales large enough to see the decay time constant, and
%> allan deviations with sampling intervals ranging from much less than the
%> time constant to sampling intervals much larger than the time constant,
%> the slope will transition from +1/2 for sampling intervals much shorter 
%> than the time constant, to 0 as the sampling interval approaches
%> the time constant at it's maxima, to -1/2 as the interval becomes much
%> larger than the time constant.
%> 
%> When multiple gauss markov processes are added together it will be
%> nearly impossible to pick apart any of these parameters, but when run
%> individually the driven noise magnitude (qc) and time constant (Tc) can be
%> found on the allan deviation plot as follows:
%> Tc = argmax_sampling_interval / 1.89, where argmax_sampling_interval is
%> the allan deviation sampling interval that maximizes the allan deviation.
%> The driven noise magnitude may then be found with the following
%> relationship:
%> qc = sigma_max / (0.437 * sqrt(Tc)), where sigma_max is the allan
%> deviation maxima.
%>
%> If the combined output of all of these Gauss-Markov processes are
%> analyzed, the only easy thing to pick off will be the velocity random
%> walk / angle random walk, which may be found by fitting a line to the
%> segment with a slope of -1/2 and finding the intersection of this line
%> with tau = 1.
%>
%> example usage:
%> gmp = GMP1('IMU1');
%> [y_out_gyro, y_out_accel, s1_gyro, b1_gyro, b2_gyro, rw_gyro,...
%>      s1_accel, b1_accel, b2_accel, rw_accel]gmp.runimu(d_theta, d_v, dt)
%>
%> Note: All values are made up but should be on the order of magnitude
%> to make sense to show 3 IMUs of differing quality.

% ======================================================================



classdef GMP1
    % GMP1 Discrete time First order Gauss-Markov process class for 
    % generating IMU errors for IMUs.  
    properties
        %> Class Public Properties
        %> Accelerometer turn on bias standard deviation
        accel_bias_stddev = ones(1,3);                % accelerometer turn on bias standard deviation
        %> Accelerometer turn on bias time constant
        accel_bias_time_const = ones(1,3);            % accelerometer turn on bias time constant (very long)
        %> Accelerometer bias inrun standard deviation
        accel_bias_inrun_stddev = ones(1,3);          % accelerometer bias inrun standard deviation
        %> Accelerometer bias inrun time constant
        accel_bias_inrun_time_const = ones(1,3);      % accelerometer bias inrun time constant (medium)
        %> Accelerometer scale factor standard deviation
        accel_scale_factor_stddev = ones(1,3);        % accelerometer scale factor standard deviation
        %> Accelerometer scale factor time constant
        accel_scale_factor_time_const = ones(1,3);    % accelerometer scale factor time constant (medium)

        %> Gyroscope turn on bias standard deviation
        gyro_bias_stddev = ones(1,3);                 % gyro turn on bias standard deviation
        %> Gyroscope turn on bias time constant
        gyro_bias_time_const = ones(1,3);             % gyro turn on bias time constant (very long)
        %> Gyroscope bias inrun standard deviation
        gyro_bias_inrun_stddev = ones(1,3);           % gyro bias inrun standard deviation
        %> Gyroscope bias inrun time constant
        gyro_bias_inrun_time_const = ones(1,3);       % gyro bias inrun time constant
        %> Gyroscope scale factor standard deviation
        gyro_scale_factor_stddev = ones(1,3);         % gyro scale factor standard deviation
        %> Gyroscope scale factor time constant
        gyro_scale_factor_time_const = ones(1,3);     % gyro scale factor time constant

        %> Gyroscope angle random walk
        gyro_angle_random_walk;                       % angle random walk
        %> Accelerometer velocity random walk
        accel_velocity_random_walk;                   % velocity random walk
        %> IMU sub-model
        imu_type;                                     % imu model number                                    
   end
   
    properties(Constant)                            % constant properties
        %> Class Constant Properties
        %> Micro g to g
        micro_g_to_g = 0.000001;
        %> Mean gravity
        gravity_mean = 9.7976446561;
        %> G to meters per second per second
        g_to_meters_per_sec2 = GMP1.gravity_mean;
        %> Hours to seconds
        hours_to_seconds = 3600;
        %> Minutes to seconds
        minutes_to_seconds = 60;
        %> Parts per million to part
        ppm_to_part = 0.000001;
        %> Degrees to radians
        deg_to_rad = pi / 180.0;
        %> Degrees per hour to radians per second
        deg_per_hr_to_rad_per_sec = GMP1.deg_to_rad / GMP1.hours_to_seconds;
        %> Per root hour to per root second
        per_root_hour_to_per_root_second = 1.0/60.0;   
   end
   
    methods
       
        % ======================================================================
        %> @brief GMP1 class constructor
        %>
        %> Instantiate a first order gauss markov process object corresponding to
        %> one of the IMU classes.  Initializes all relevant properties to
        %> those specified by the <redacted> source file.
        %>
        %> @param arg Sensor model, currently accepted inputs are 'IMU1',
        %> 'IMU2', and 'IMU3'.
        %>
        %> @return instance of the GMP1 class
        % ======================================================================
        function this = GMP1(arg)         
              % Class constructor
              import constants.*;
              if strcmp(arg, 'IMU1')
                  this.imu_type = 'IMU1';
                  fprintf('IMU configured to %s \n', arg);

                  % accelerometer coefficients
                  %> Set the accelerometer turn on bias standard deviation
                  %> convert from micro g to m/s/s.
                  this.accel_bias_stddev(1) = 10800.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_stddev(2) = 10800.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_stddev(3) = 10800.0*this.micro_g_to_g * this.g_to_meters_per_sec2;

                  %> Set the accelerometer turn on bias standard deviation
                  %> convert from hours to seconds
                  this.accel_bias_time_const(1) = 1080.0*this.hours_to_seconds;
                  this.accel_bias_time_const(2) = 1080.0*this.hours_to_seconds;
                  this.accel_bias_time_const(3) = 1080.0*this.hours_to_seconds;

                  %> Set the accelerometer bias inrun standard deviation
                  %> convert from micro g to m/s/s
                  this.accel_bias_inrun_stddev(1) = 400.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_inrun_stddev(2) = 400.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_inrun_stddev(3) = 400.0*this.micro_g_to_g * this.g_to_meters_per_sec2;

                  %> Set the accelerometer bias inrun time constant
                  %> convert from minutes to seconds
                  this.accel_bias_inrun_time_const(1) = 5.5 * this.minutes_to_seconds;
                  this.accel_bias_inrun_time_const(2) = 5.5 * this.minutes_to_seconds;
                  this.accel_bias_inrun_time_const(3) = 5.5 * this.minutes_to_seconds;

                  %> Set the accelerometer scale factor standard deviation
                  %> convert from parts per million to parts
                  this.accel_scale_factor_stddev(1) = 1110 * this.ppm_to_part;
                  this.accel_scale_factor_stddev(2) = 1110 * this.ppm_to_part;
                  this.accel_scale_factor_stddev(3) = 1110 * this.ppm_to_part;

                  %> Set the accelerometer scale factor time constant
                  %> convert from minutes to seconds
                  this.accel_scale_factor_time_const(1) = 5.5 * this.minutes_to_seconds;
                  this.accel_scale_factor_time_const(2) = 5.5 * this.minutes_to_seconds;
                  this.accel_scale_factor_time_const(3) = 5.5 * this.minutes_to_seconds;

                  %> set the accelerometer velocity random walk standard dev
                  %> convert from m/s/sqrt(hr) to m/s/sqrt(s)
                  this.accel_velocity_random_walk = 0.152 * this.per_root_hour_to_per_root_second;

    %               accel_SF_nonlinearity.r(0)->Set_initial_stdv( 178.0 * CONSTANTS::ppm_to_part/CONSTANTS::gravity_mean ); //ppm/g to parts/(m/s^2)
    %               accel_SF_nonlinearity.r(0)->Set_time_constant( 5.5 * CONSTANTS::minutes_to_seconds );
    % 
    %               accel_nonorthogonality.r(0)->Set_initial_stdv( 170.0 * CONSTANTS::arcseconds_to_radians );
    %               accel_nonorthogonality.r(0)->Set_time_constant( 5.5 * CONSTANTS::minutes_to_seconds );              
    %               
    %               accel_misalign.r(0)->Set_initial_stdv( 170.0 * CONSTANTS::arcseconds_to_radians );
    %               accel_misalign.r(0)->Set_time_constant( 5.5 * CONSTANTS::minutes_to_seconds );

                  % gyro coefficients
                  %> Set the gyro turn on bias standard deviation
                  %> convert from deg per hour to rad per sec
                  this.gyro_bias_stddev(1) = 70.0 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_stddev(2) = 70.0 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_stddev(3) = 70.0 * this.deg_per_hr_to_rad_per_sec;

                  %> Set the gyro turn on bias time constant
                  %> convert from hours to seconds
                  this.gyro_bias_time_const(1) = 1100.0 * this.hours_to_seconds;
                  this.gyro_bias_time_const(2) = 1100.0 * this.hours_to_seconds;
                  this.gyro_bias_time_const(3) = 1100.0 * this.hours_to_seconds;

                  %> Set the gyro bias in run standard deviation
                  %> convert from deg per hour to rad per sec
                  this.gyro_bias_inrun_stddev(1) = 1.4 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_inrun_stddev(2) = 1.4 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_inrun_stddev(3) = 1.4 * this.deg_per_hr_to_rad_per_sec;

                  %> Set the gyro bias inrun time constant
                  %> convert from minutes to seconds
                  this.gyro_bias_inrun_time_const(1) = 5.4 * this.minutes_to_seconds;
                  this.gyro_bias_inrun_time_const(2) = 5.4 * this.minutes_to_seconds;
                  this.gyro_bias_inrun_time_const(3) = 5.4 * this.minutes_to_seconds;

                  %> Set the gyro scale factor standard deviation
                  %> convert from parts per million to parts
                  this.gyro_scale_factor_stddev(1) = 1120.0 * this.ppm_to_part;
                  this.gyro_scale_factor_stddev(2) = 1120.0 * this.ppm_to_part;
                  this.gyro_scale_factor_stddev(3) = 1120.0 * this.ppm_to_part;

                  %> Set the gyro scale factor time constant
                  %> convert from minutes to seconds
                  this.gyro_scale_factor_time_const(1) = 5.4 * this.minutes_to_seconds;
                  this.gyro_scale_factor_time_const(2) = 5.4 * this.minutes_to_seconds;
                  this.gyro_scale_factor_time_const(3) = 5.4 * this.minutes_to_seconds;

                  %> Set the gyro angle random walk standard deviation
                  %> convert from deg per root hour to rand per root second
                  this.gyro_angle_random_walk = 0.178 * this.deg_to_rad * this.per_root_hour_to_per_root_second;


    %               gyro_nonorthogonality.r(0)->Set_initial_stdv( 185.0* CONSTANTS::arcseconds_to_radians );
    %               gyro_nonorthogonality.r(0)->Set_time_constant( 5.4 * CONSTANTS::minutes_to_seconds );


              elseif strcmp(arg, 'IMU2')
                  this.imu_type = 'IMU2';
                  fprintf('IMU configured to %s \n', arg);

                  % accel
                  this.accel_bias_stddev(1) = 10800.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_stddev(2) = 10800.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_stddev(3) = 10800.0*this.micro_g_to_g * this.g_to_meters_per_sec2;

                  this.accel_bias_time_const(1) = 1080.0*this.hours_to_seconds;
                  this.accel_bias_time_const(2) = 1080.0*this.hours_to_seconds;
                  this.accel_bias_time_const(3) = 1080.0*this.hours_to_seconds;

                  this.accel_bias_inrun_stddev(1) = 400.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_inrun_stddev(2) = 400.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_inrun_stddev(3) = 400.0*this.micro_g_to_g * this.g_to_meters_per_sec2;

                  this.accel_bias_inrun_time_const(1) = 5.5 * this.minutes_to_seconds;
                  this.accel_bias_inrun_time_const(2) = 5.5 * this.minutes_to_seconds;
                  this.accel_bias_inrun_time_const(3) = 5.5 * this.minutes_to_seconds;

                  this.accel_scale_factor_stddev(1) = 1110 * this.ppm_to_part;
                  this.accel_scale_factor_stddev(2) = 1110 * this.ppm_to_part;
                  this.accel_scale_factor_stddev(3) = 1110 * this.ppm_to_part;

                  this.accel_scale_factor_time_const(1) = 5.5 * this.minutes_to_seconds;
                  this.accel_scale_factor_time_const(2) = 5.5 * this.minutes_to_seconds;
                  this.accel_scale_factor_time_const(3) = 5.5 * this.minutes_to_seconds;

                  this.accel_velocity_random_walk = 0.089 * this.per_root_hour_to_per_root_second;

    %               accel_SF_nonlinearity.r(0)->Set_initial_stdv( 152.0 * CONSTANTS::ppm_to_part/CONSTANTS::gravity_mean ); //ppm/g to parts/(m/s^2)
    %               accel_SF_nonlinearity.r(0)->Set_time_constant( 5.4 * CONSTANTS::minutes_to_seconds );
    %          
    %               accel_nonorthogonality.r(0)->Set_initial_stdv( 152.0 * CONSTANTS::arcseconds_to_radians );
    %               accel_nonorthogonality.r(0)->Set_time_constant( 5.4 * CONSTANTS::minutes_to_seconds );
    %          
    %               accel_misalign.r(0)->Set_initial_stdv( 152.0 * CONSTANTS::arcseconds_to_radians );
    %               accel_misalign.r(0)->Set_time_constant( 5.4 * CONSTANTS::minutes_to_seconds );

                  % gyro
                  this.gyro_bias_stddev(1) = 39.2 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_stddev(2) = 39.2 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_stddev(3) = 39.2 * this.deg_per_hr_to_rad_per_sec;

                  this.gyro_bias_time_const(1) = 1100.0 * this.hours_to_seconds;
                  this.gyro_bias_time_const(2) = 1100.0 * this.hours_to_seconds;
                  this.gyro_bias_time_const(3) = 1100.0 * this.hours_to_seconds;

                  this.gyro_bias_inrun_stddev(1) = 1.3 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_inrun_stddev(2) = 1.3 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_inrun_stddev(3) = 1.3 * this.deg_per_hr_to_rad_per_sec;

                  this.gyro_bias_inrun_time_const(1) = 5.2 * this.minutes_to_seconds;
                  this.gyro_bias_inrun_time_const(2) = 5.2 * this.minutes_to_seconds;
                  this.gyro_bias_inrun_time_const(3) = 5.2 * this.minutes_to_seconds;

                  this.gyro_scale_factor_stddev(1) = 912.0 * this.ppm_to_part;
                  this.gyro_scale_factor_stddev(2) = 912.0 * this.ppm_to_part;
                  this.gyro_scale_factor_stddev(3) = 912.0 * this.ppm_to_part;

                  this.gyro_scale_factor_time_const(1) = 5.3 * this.minutes_to_seconds;
                  this.gyro_scale_factor_time_const(2) = 5.3 * this.minutes_to_seconds;
                  this.gyro_scale_factor_time_const(3) = 5.3 * this.minutes_to_seconds;

                  this.gyro_angle_random_walk = 0.135 * this.deg_to_rad * this.per_root_hour_to_per_root_second;

    %               gyro_nonorthogonality.r(0)->Set_initial_stdv( 142.0* CONSTANTS::arcseconds_to_radians );
    %               gyro_nonorthogonality.r(0)->Set_time_constant( 5.4 * CONSTANTS::minutes_to_seconds );

              elseif strcmp(arg, 'IMU3')
                  this.imu_type = 'IMU3';
                  fprintf('IMU configured to %s \n', arg);

                  % accel
                  this.accel_bias_stddev(1) = 5310.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_stddev(2) = 5310.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_stddev(3) = 5310.0*this.micro_g_to_g * this.g_to_meters_per_sec2;

                  this.accel_bias_time_const(1) = 991.0*this.hours_to_seconds;
                  this.accel_bias_time_const(2) = 991.0*this.hours_to_seconds;
                  this.accel_bias_time_const(3) = 991.0*this.hours_to_seconds;

                  this.accel_bias_inrun_stddev(1) = 270.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_inrun_stddev(2) = 270.0*this.micro_g_to_g * this.g_to_meters_per_sec2;
                  this.accel_bias_inrun_stddev(3) = 270.0*this.micro_g_to_g * this.g_to_meters_per_sec2;

                  this.accel_bias_inrun_time_const(1) = 5.4 * this.minutes_to_seconds;
                  this.accel_bias_inrun_time_const(2) = 5.4 * this.minutes_to_seconds;
                  this.accel_bias_inrun_time_const(3) = 5.4 * this.minutes_to_seconds;

                  this.accel_scale_factor_stddev(1) = 740.0 * this.ppm_to_part;
                  this.accel_scale_factor_stddev(2) = 740.0 * this.ppm_to_part;
                  this.accel_scale_factor_stddev(3) = 740.0 * this.ppm_to_part;

                  this.accel_scale_factor_time_const(1) = 5.4 * this.minutes_to_seconds;
                  this.accel_scale_factor_time_const(2) = 5.4 * this.minutes_to_seconds;
                  this.accel_scale_factor_time_const(3) = 5.4 * this.minutes_to_seconds;

                  this.accel_velocity_random_walk = 0.089 * this.per_root_hour_to_per_root_second;

    %               accel_SF_nonlinearity.r(0)->Set_initial_stdv( 151.0 * CONSTANTS::ppm_to_part/CONSTANTS::gravity_mean ); //ppm/g to parts/(m/s^2)
    %               accel_SF_nonlinearity.r(0)->Set_time_constant( 5.2 * CONSTANTS::minutes_to_seconds );
    %          
    %               accel_nonorthogonality.r(0)->Set_initial_stdv( 101.1 * CONSTANTS::arcseconds_to_radians );
    %               accel_nonorthogonality.r(0)->Set_time_constant( 5.2 * CONSTANTS::minutes_to_seconds );
    %          
    %               accel_misalign.r(0)->Set_initial_stdv( 151.0 * CONSTANTS::arcseconds_to_radians );
    %               accel_misalign.r(0)->Set_time_constant( 5.2 * CONSTANTS::minutes_to_seconds );

                  % gyro
                  this.gyro_bias_stddev(1) = 17.5 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_stddev(2) = 17.5 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_stddev(3) = 17.5 * this.deg_per_hr_to_rad_per_sec;

                  this.gyro_bias_time_const(1) = 920.0 * this.hours_to_seconds;
                  this.gyro_bias_time_const(2) = 920.0 * this.hours_to_seconds;
                  this.gyro_bias_time_const(3) = 920.0 * this.hours_to_seconds;

                  this.gyro_bias_inrun_stddev(1) = 1.1 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_inrun_stddev(2) = 1.1 * this.deg_per_hr_to_rad_per_sec;
                  this.gyro_bias_inrun_stddev(3) = 1.1 * this.deg_per_hr_to_rad_per_sec;

                  this.gyro_bias_inrun_time_const(1) = 5.5 * this.minutes_to_seconds;
                  this.gyro_bias_inrun_time_const(2) = 5.5 * this.minutes_to_seconds;
                  this.gyro_bias_inrun_time_const(3) = 5.5 * this.minutes_to_seconds;

                  this.gyro_scale_factor_stddev(1) = 600.0 * this.ppm_to_part;
                  this.gyro_scale_factor_stddev(2) = 600.0 * this.ppm_to_part;
                  this.gyro_scale_factor_stddev(3) = 600.0 * this.ppm_to_part;

                  this.gyro_scale_factor_time_const(1) = 5.5 * this.minutes_to_seconds;
                  this.gyro_scale_factor_time_const(2) = 5.5 * this.minutes_to_seconds;
                  this.gyro_scale_factor_time_const(3) = 5.5 * this.minutes_to_seconds;

                  this.gyro_angle_random_walk = 0.120 * this.deg_to_rad * this.per_root_hour_to_per_root_second;

    %               gyro_nonorthogonality.r(0)->Set_initial_stdv( 100.0* CONSTANTS::arcseconds_to_radians );
    %               gyro_nonorthogonality.r(0)->Set_time_constant( 5.5 * CONSTANTS::minutes_to_seconds );

              else
                  error('Error: unsupported imu model, must be IMU1, IMU2, or IMU3');
              end
        end

        % ======================================================================
        %> @brief run_gauss_markov generate first order gauss markov time series 
        %>
        %> Generate a first order Gauss-Markov time series for scale factor,
        %> turn on bias, and bias inrun, as well as a random walk time series.
        %> Returns the individual processes, as well as the combined output to
        %> lay errors on top of synthetic, perfect, data.
        %>
        %> @param this instance of the GMP1 class.
        %> @param sensor_type sensor_type that the user would like to generate
        %> time series for.  Options are 0 = gyroscope, 1 = accelerometer.
        %> @param y_true time series that errors are to be generated on top of.
        %> @param dt sampling time of the process
        %>
        %> @retval y_out output time series vector with overlaid simulated noise
        %> @retval s1 scale factor time series vector (for debug / validation)
        %> @retval b1 bias inrun time series vector (for debug / validation)
        %> @retval b2 turn on bias time series vector (for debug / validation)
        %> @retval rw random walk time series vector (for debug / validation)
        % ======================================================================

        function [y_out, s1, b1, b2, rw] = run_gauss_markov(this, sensor_type, y_true, dt)
                tic;

                len = length(y_true);

                % scale factor gauss-markov process vector
                s1 = ones(len,1);

                % bias term 1 gauss-markov process vector
                b1 = ones(len,1);

                % bias term 2 gauss-markov process vector
                b2 = ones(len,1);

                mu = 0; % Gaussian noise
                random_walk_vector = ones(len,1);

                % set time constant terms
                % note: 1/tau = correlation time

                if sensor_type == 0
                    fprintf('Processing gyroscope\n')
                    % set time constant terms
                    tau_s1 = this.gyro_scale_factor_time_const(1);      % sec
                    tau_b1 = this.gyro_bias_inrun_time_const(1);        % sec
                    tau_b2 = this.gyro_bias_time_const(1);              % sec

                    % set standard deviation terms
                    sigma_s1 = this.gyro_scale_factor_stddev(1);        % unitless
                    sigma_b1 = this.gyro_bias_inrun_stddev(1);          % rad/s
                    sigma_b2 = this.gyro_bias_stddev(1);                % rad/s

                    random_walk_sigma = this.gyro_angle_random_walk;          % rad/rt-sec

                    for i=1:length(random_walk_vector)
                       random_walk_vector(i) = normrnd(mu, random_walk_sigma); 
                    end            

                elseif sensor_type == 1
                    fprintf('Processing accelerometer\n')
                    % set time constant terms
                    tau_s1 = this.accel_scale_factor_time_const(1);     % sec
                    tau_b1 = this.accel_bias_inrun_time_const(1);       % sec
                    tau_b2 = this.accel_bias_time_const(1);             % sec

                    % set standard deviation terms
                    sigma_s1 = this.accel_scale_factor_stddev(1);       % m/s/s
                    sigma_b1 = this.accel_bias_inrun_stddev(1);         % m/s/s
                    sigma_b2 = this.accel_bias_stddev(1);               % m/s/s

                    random_walk_sigma = this.accel_velocity_random_walk;      % m/s/rt-sec

                    for i=1:length(random_walk_vector)
                       random_walk_vector(i) = normrnd(mu, random_walk_sigma); 
                    end

                else
                    fprintf('Error: invalid sensor type.\n');
                end

                % initialize markov process kernel to a gaussian sample from
                % the specified standard deviations

                s_init = normrnd(mu, sigma_s1);
                b1_init = normrnd(mu, sigma_b1);
                b2_init = normrnd(mu, sigma_b2);

                % scale initial values by dt
                s1(1) = s_init * dt;
                b1(1) = b1_init * dt;
                b2(1) = b2_init * dt;

                % run the gauss markov process for scalefactor and bias terms
                for i = 2:len 
                        s1(i) = (exp((-1/tau_s1) * dt) * s1(i-1) + normrnd(mu,sigma_s1) * dt); % defined by sf time const and std dev
                        b1(i) = (exp((-1/tau_b1) * dt) * b1(i-1) + normrnd(mu,sigma_b1) * dt); % defined by bias inrun time const and std dev  
                        b2(i) = (exp((-1/tau_b2) * dt) * b2(i-1) + normrnd(mu,sigma_b2) * dt); % defined by bias time const and std dev  
                end

                % note scaling bias terms for dimensions to work
                % b1 and b2 are in m/s/s, so multiplication by dt gives delta
                % m/s, and VRW is in m/s/sqrt(s), so scaling by sqrt(s)
                % achieves dimensional consistency

                scale = s1 + 0; % no gaussian term for scale factor
                b1 = b1.*dt;            % inrun
                b2 = b2.*dt;            % quasi static (turn on)
                rw = random_walk_vector.*(dt^0.5);

                bias = b1 + b2 + rw;

                y_out = y_true .* (1 + scale) + bias;

                toc;
        end

        % ======================================================================
        %> @brief run_imu method to generate 1-D IMU errors
        %>
        %> @param this instance of the GMP1 class
        %> @param y_true_gyro error-free 1D gyro input vector (nx1)
        %> @param y_true_accel error-free 1D accelerometer input vector (nx1)
        %> @param dt sampling time (seconds)
        %>
        %> @ret returns concatenation of output vectors from run_gauss_markov
        % ======================================================================

        function [y_out_gyro, y_out_accel, s1_gyro, b1_gyro, b2_gyro, rw_gyro,...
                    s1_accel, b1_accel, b2_accel, rw_accel] = run_imu(this, y_true_gyro, y_true_accel, dt)

                [y_out_gyro, s1_gyro, b1_gyro, b2_gyro, rw_gyro] = this.run_gauss_markov(0, y_true_gyro, dt);
                [y_out_accel, s1_accel, b1_accel, b2_accel, rw_accel] = this.run_gauss_markov(1, y_true_accel, dt);          
        end


        % ======================================================================
        %> @brief run_imu3d method to generate 3-D IMU errors
        %>
        %> @param this instance of the GMP1 class
        %> @param y_true_gyro error-free 3D gyro input matrix (nx3)
        %> @param y_true_accel error-free 3D accelerometer input matrix (nx3)
        %> @param dt sampling time (seconds)
        %>
        %> @ret returns concatenation of output vectors from run_gauss_markov
        %> concatenated into nx3 matrices
        % ======================================================================
        function [y_out_gyro, y_out_accel, s1_gyro, b1_gyro, b2_gyro, rw_gyro,...
                    s1_accel, b1_accel, b2_accel, rw_accel] = run_imu3d(this, y_true_gyro, y_true_accel, dt)

            for i = 1:3
                [y_out_gyro{i}, s1_gyro{i}, b1_gyro{i}, b2_gyro{i}, rw_gyro{i}] = this.run_gauss_markov(0, y_true_gyro(:,i), dt);
                [y_out_accel{i}, s1_accel{i}, b1_accel{i}, b2_accel{i}, rw_accel{i}] = this.run_gauss_markov(1, y_true_accel(:,i), dt);
            end
            y_out_gyro = cell2mat(y_out_gyro);
            y_out_accel = cell2mat(y_out_accel);

            s1_gyro = cell2mat(s1_gyro);
            b1_gyro = cell2mat(b1_gyro);
            b2_gyro = cell2mat(b2_gyro);
            rw_gyro = cell2mat(rw_gyro);

            s1_accel = cell2mat(s1_accel);
            b1_accel = cell2mat(b1_accel);
            b2_accel = cell2mat(b2_accel);
            rw_accel = cell2mat(rw_accel);

        end


        % ======================================================================
        %> @brief print_coeffs method to print out Gauss-Markov process
        %> coefficients
        %>
        %> Used to verify coefficients, printed in units consistent with ECTOS,
        %> not consistent with normal human readable (datasheet) units.
        % ======================================================================
        function print_coeffs(this)
            fprintf('Gyro Scale Factor Time Constant: %f\n', this.gyro_scale_factor_time_const(1));
            fprintf('Gyro Bias Inrun Time Constant: %f\n', this.gyro_bias_inrun_time_const(1));
            fprintf('Gyro Bias Time Constant: %f\n', this.gyro_bias_time_const(1));
            fprintf('Gyro Scale Factor stddev: %f\n', this.gyro_scale_factor_stddev(1));
            fprintf('Gyro Bias Inrun stddev: %f\n', this.gyro_bias_inrun_stddev(1));
            fprintf('Gyro Bias stddev: %f\n', this.gyro_bias_stddev(1));
            fprintf('Gyro Angle Random Walk: %f\n\n', this.gyro_angle_random_walk);

            fprintf('Accel Scale Factor Time Constant: %f\n', this.accel_scale_factor_time_const(1));
            fprintf('Accel Bias Inrun Time Constant: %f\n', this.accel_bias_inrun_time_const(1));
            fprintf('Accel Bias Time Constant: %f\n', this.accel_bias_time_const(1));
            fprintf('Accel Scale Factor stddev: %f\n', this.accel_scale_factor_stddev(1));
            fprintf('Accel Bias Inrun stddev: %f\n', this.accel_bias_inrun_stddev(1));
            fprintf('Accel Bias stddev: %f\n', this.accel_bias_stddev(1));
            fprintf('Accel Velocity Random Walk: %f\n\n', this.accel_velocity_random_walk);
        end
      
      
   % end methods  
   end
% end class
end
