%> @file BaroErrors.m
%> @brief Class for adding barometeric altimeter errors to time series data
%>
%> @author Eric Victorson
%> @date 2018/10/15
%> @version 1.0
% ======================================================================
%> @brief BaroErrors class to model barometric altimeter errors in time
%> series data
%>
%> Currently only modelled as a uniformly distributed turn on bias that
%> defaults to -x to x meters if no input is supplied, with an
%> overlying white noise measurement process that defaults to a standard
%> deviation of x m.
% ======================================================================

classdef BaroErrors
 
    %> Class public properties.  Only white noise supported for now.
    properties
        %> standard deviation for measurement white noise process
        sigma;
        %> uniformly distributed turn on bias for altimeter
        bias;
    end
    
    
    methods
        
        % ======================================================================
        %> @brief BaroErrors class constructor
        %>
        %> Instantiate a BaroErrors object to add white noise and bias to
        %> barometric altimeter measurements. 0.1 m seems reasonable.
        %>
        %> @param varargin{1} standard deviation for white noise
        %> @param varargin{2} uniform distribution maxima for bias
        %>
        %> @return this instance of the BaroErrors class
        % ======================================================================
        function this = BaroErrors(varargin)
            if nargin > 0
                this.sigma = varargin{1};
                this.bias = (-1 + rand(1)*2)*varargin;
            else
								bias_mag = 10;
                this.sigma = 0.1;
                this.bias = rand(120);
                this.bias = (-1 + rand(1)*2)*bias_mag;
            end
        end
        
        % ======================================================================
        %> @brief add_bias add sample from uniform distribution as bias
        %>
        %> @param yin error free altimeter time series
        %>
        %> @retval y altimeter time series with added bias
        % ======================================================================
        function y = add_bias(this, yin)
            y = yin + this.bias;
        end
        
        % ======================================================================
        %> @brief add_white_noise add white noise to time series
        %>
        %> @param yin input time series
        %>
        %> @retval y time series with added white noise 
        % ======================================================================
        function y = add_white_noise(this, yin)
           len = length(yin);
           y = zeros(len, 1);
           for i = 1:len
              y(i) = yin(i) + normrnd(0, this.sigma); 
           end
           
        end
        
        % ======================================================================
        %> @brief process_altimeter add bias and white noise to altimeter
        %> data
        %>
        %> @param yin input time series
        %>
        %> @retval y output time series with added bias and white noise
        % ======================================================================
        function y = process_altimeter(this, yin)
           y_temp = this.add_bias(yin);
           y = this.add_white_noise(y_temp);
        end
        
        % end methods
    end
    
    % end class
end
