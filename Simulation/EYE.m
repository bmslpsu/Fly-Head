classdef EYE
    % EYE: represents the spatiotemporal attributes of an animals visual
    % system
    
    properties
        delta_phi       % angle between adjacent ommatidia [deg]
        time_constant   % temporal filter time constant [s]
        delta_rho       % 
        filt            % spatial blurring filter
       	n_receptor     	% # of ommatidia (default = 2)
        n_pts           % sample points
    end
    
    methods
        function obj = EYE(delta_phi,time_constant,n)
            % EYE Construct an instance of this class
            %   Input the angle between adjacent ommatidia & the # of ommatidia
            if nargin==0
                delta_phi = 4.6*pi/180; % default for Drosophila (rough approximation; follows from 
                % caption of Fig. 18, Buchner, 1981 (in Ali))
                n = 2;
            elseif nargin==1
                n = 2;
            end
            
            obj.delta_phi = delta_phi; % angle between adjacent ommatidia
            obj.time_constant = time_constant;
            obj = SetProp(obj,n);
       	end
        
        function obj = SetProp(obj,n)
            % SetProp: set properties of EYE
            %   Use the angle between adjacent ommatidia % # of ommatidia
            %   to compute properties of EYE
            
          	obj.n_receptor  = n;                        % # of ommatidia
            obj.delta_rho   = obj.delta_phi*1.1;        % 
            theta           = -pi:pi/480:pi - pi/480;  	% 

            % From Snyder (1979) as cited in Burton & Laughlin (2003)
            obj.filt =  exp( -4.*log(2).*abs(theta).^2 ./ obj.delta_rho^2 ); % spatial blurring filter
            
            mid = obj.n_receptor/2 ;
            
            eye_filt(:,mid+1) = obj.filt;
            cnt = 1;
            for j = (mid+2):obj.n_receptor
                eye_filt(:,j) = circshift(obj.filt, [0 cnt*15]);
                cnt = cnt + 1;
            end
            cnt = 1;
            for j = mid:-1:1
                eye_filt(:,j) = circshift(obj.filt, [0 -cnt*15]);
                cnt = cnt + 1;
            end
            
            obj.filt = eye_filt;

            eye_filt(eye_filt < 0.005) = 0;  % set vey low values to
            % zero, so may be able to use eye_filt as a sparse matrix...not a time
            % saver, because gets mulitples with non-sparse matrices.

            [obj.n_pts, obj.n_receptor] = size(eye_filt);
        end
    end
end

