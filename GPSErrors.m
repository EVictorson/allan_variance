%> @file GPSErrors.m
%> @brief Class for adding GPS errors to time series
%>
%> @author Eric Victorson
%> @date 2018/10/15
%> @version 1.0
% ======================================================================
%> @brief GPSErrors class to model gps errors in time series data
%>
%> Neglect atmospheric GPS errors and just use simple gaussian white noise
%> for now.  ~1-3m standard deviation is reasonable.
%>
%> Contributions to GPS errors include:
%>
%> Ionospheric effects: clouds of free elections that act as a dispersive 
%>  medium for GPS signals.  The primary effect is to change the signal
%>  propogation speed as compared to that of free space.  The signal modulation
%>  is delayed while the carrier phase is advanced by the same amount.  
%> Ephemeris errors: Small errors in the ephemeris data transmitted by each
%   satellite
%> Satellite clock errors: Each satellite atomic clock is allowed some
%   degree of relative drift that is estimated by a ground station, which is
%   used to generate clock correction data in the GPS navigation message.
%   When space vehicle (SV) time is corrected using this data the result is
%   called GPS time.
%> Multipath distortion: Signal reflections off of objects causing a
%   lengthening of the propogation path
%> Tropospheric effects: lengthens the propogation path due to
%>  refractions with dry gases and water vapor.
%>
%> All of the above errors can be converted to into an equivalent range
%> error experienced by the user, called the user-equivalent range error
%> (UERE)
% ======================================================================


classdef GPSErrors
 
    %> Class public properties.  Only white noise supported for now.
    properties
        %> standard deviation for measurement white noise process
        sigma;
        
    end
    
    
    methods
        
        % ======================================================================
        %> @brief GPSErrors class constructor
        %>
        %> Instantiate a GPSErrors object to add white noise to GPS
        %> measurements
        %> @param varargin standard deviation for white noise
        %>
        %> @return this instance of GPSErrors class
        % ======================================================================
        function this = GPSErrors(varargin)
            if nargin > 0 
                this.sigma = varargin{1};
            else
                this.sigma = 3;
            end
        end

        % ======================================================================
        %> @brief add_white_noise function to add white noise to GPS time
        %> series
        %>
        %> @param y_true error free gps time series
        %>
        %> @retval y gps time series with gaussian noise
        % ======================================================================
        function y = add_white_noise(this, y_true)
           len = length(y_true);
           y = zeros(len,1);
           for i = 1:len
              y(i) = y_true(i) + normrnd(0, this.sigma); 
           end
            
        end
        
        % ======================================================================
        %> @brief process_gps3d add white noise to 3d gps matrix
        %>
        %> @param y_true3d error free 3d gps matrix
        %>
        %> @retval y3d 3d gps time series with gaussian noise
        % ======================================================================
        function y3d = process_gps3d(this, y_true3d)
            for i = 1:3
                y3d_temp{i} = this.add_white_noise(y_true3d(:,i)); 
            end
            y3d = cell2mat(y3d_temp);
        end
        
        
       % end methods 
    end
    % end class
end
