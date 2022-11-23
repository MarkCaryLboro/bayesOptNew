classdef aei < acqFcn
    % Implement the eadaptive xpected improvement acquisition function

    properties ( Constant = true )
        FcnName         acqFcnType   = acqFcnType( "aei" )
    end % Abstract constant properties
    
    properties ( SetAccess = protected )
        Beta  
        MnExpRate  (1,:)    double     = 0.25
        MxExpRate  (1,:)    double     = 2.00
    end

    properties ( Access = protected )
        Count      (1,1)    int64     { mustBeGreaterThanOrEqual( Count,... % Counter for determining when to adjust memory rate
                                                                  0 ) } = 0
        Reset      (1,1)    logical = true
    end % protected properties

    properties ( SetAccess = protected, Dependent = true )
        R                                                                   % Exploration rate
    end % Dependent properties

    methods
        function obj = aei( ModelObj, Beta )
            %--------------------------------------------------------------
            % Adaptive Expected Improvement (AEI) class constructor
            %
            % obj = aei( ModelObj, Beta );
            %
            % Input Arguments:
            %
            % ModelObj    --> Surrogate model for system
            % Beta        --> (double) Hyperparameter {0.01}
            %--------------------------------------------------------------
            arguments
                ModelObj (1,1)         { mustBeNonempty( ModelObj ) }
                Beta     (1,3)  double      = [ 1.0, 1.05, 0.05 ]           % [ R0, M, D ]
            end
            obj = obj.setBeta( Beta );
            obj.ModelObj = ModelObj;
        end % constructor

        function obj = setExpRateLimits( obj, MnRate, MxRate )
            %--------------------------------------------------------------
            % Set the minimum and maximum exploration rates
            %
            % obj = obj.setExpRateLimits( MnRate, MxRate );
            %
            % Input Arguments:
            %
            % MnRate    --> (double) Minimum exploration rate
            % MxRate    --> (double) Maximum exploration rate
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   aei     { mustBeNonempty( obj ) }
                MnRate  (1,1)   double  { mustBeGreaterThan( MnRate, 0 ) } = obj.MnRate
                MxRate  (1,1)   double  { mustBeLessThan( MxRate, 2.5 ) }  = obj.MxRate
            end
            try                                                             %#ok<TRYNC> 
                obj.MnExpRate = MnRate;
            end
            try                                                             %#ok<TRYNC> 
                obj.MxExpRate = MxRate;
            end
        end % setExpRateLimits

        function Fcn = evalFcn( obj, X, Beta )
            %--------------------------------------------------------------
            % Evaluate the EI acquisition function at the location
            % specified
            %
            % Fcn = obj.evalFcn( X, Beta );
            %
            % Input Arguments:
            %
            % X       --> Points to evaluate the acquisition function
            % Beta    --> Hyperparameter value
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   aei         {mustBeNonempty( obj )}
                X       (:,:)   double      {mustBeNonempty( X )}
                Beta    (1,3) double      = obj.Beta
            end
            %--------------------------------------------------------------
            % Need a persitent variable for the rate update
            %--------------------------------------------------------------
            persistent Xlast
            if isempty( Xlast ) || obj.Reset
                Xlast = X;
                obj.Reset = false;
            end
            obj = obj.setBeta( Beta );
            [ ~, ~, D ] = obj.getHyperParams();
            Xc = obj.ModelObj.code( X );
            Xstar = obj.ModelObj.code( obj.BestX );
            CodedXlast = obj.ModelObj.code( Xlast );
            if ( norm( Xc - Xstar, 2 ) < D ) &&...
               ( norm( Xc - CodedXlast ) > 1e-4 )
                %----------------------------------------------------------
                % Increase the exploration rate
                %----------------------------------------------------------
                obj.Count = obj.Count + 1;
            end
            Xlast = X;
            %--------------------------------------------------------------
            % Calculate the acquisition function
            %--------------------------------------------------------------
            [ Z, Mu, Sigma ] = obj.calcZscore( X );
            Zpdf = normpdf( Z );
            Zcdf = normcdf( Z );
            Fcn = ( Mu - obj.Fmax ) .* Zcdf + obj.R * Sigma .* Zpdf;
            Fcn = - Fcn;
            Fcn = ( Sigma > 0 ) .* Fcn;
        end % evalFcn

        function obj = setBeta( obj, Beta )
            %--------------------------------------------------------------
            % Set the hyperparameter
            %
            % obj = obj.setBeta( Beta )
            %
            % Input Arguments:
            %
            % Beta    --> (double) hyper-parameter vector
            %--------------------------------------------------------------
            arguments
                obj  (1,1)      aei       {mustBeNonempty( obj )}
                Beta (1,3)      double  = [ 1, 1.05, 0.25 ]                 % [ R0, M, D ]
            end
            Ok = ( Beta >= [ 0.1, 1, sqrt(eps) ] );
            assert( all( Ok ), "Hyperparameter outside minimium bound" );
            obj.Beta = Beta;
        end % setBeta

        function [ Z, Mu, Sigma ] = calcZscore( obj, X )
            %--------------------------------------------------------------
            % Calculate the adjusted standard normal variable
            %
            % [ Z, Mu, Sigma ] = obj.calcZscore( X );
            %
            % Input Arguments:
            %
            % X         --> (double) Input matrix
            %
            % Output Arguments:
            %
            % Z     - Z-scores
            % Mu    - Surrogate model evaluations g(X).
            % Sigma - Surrogate model standard deviations
            %--------------------------------------------------------------
            arguments
                obj (1,1)   aei         { mustBeNonempty( obj ) }
                X   (:,:)   double      { mustBeNonempty( X ) }                
            end
            [ Mu, Sigma ] = obj.ModelObj.predict( X );
            try
                Z = ( Mu - obj.ModelObj.Fmax ) ./ Sigma;
            catch 
                Z = zeros( size( Mu ) );
            end
        end % calcZscore

        function obj = resetCounter( obj )
            %--------------------------------------------------------------
            % Reset the exploration rate index counter to zero
            %
            % obj = obj.resetCounter();
            %--------------------------------------------------------------
            obj.Count = 0;
            obj.Reset = true;
        end % resetCounter
    end % ordinary method signatures

    methods
        function R = get.R( obj )
            % Return the exploration rate
            [ R0, M ] = obj.getHyperParams();
            R = R0 * M .^ double( obj.Count );
            %--------------------------------------------------------------
            % Apply clips
            %--------------------------------------------------------------
            R = max( [ R, obj.MnExpRate] );
            R = min( [ R, obj.MxExpRate ] );
        end % get.R
    end % Set/Get methods

    methods ( Access = private )
        function [ R0, M, D] = getHyperParams( obj )
            %--------------------------------------------------------------
            % Return the hyperparameters
            %
            % [ R0, M, D] = obj.getHyperParams();
            %
            % Output Arguments:
            %
            % R0    --> Initial exploration rate
            % M     --> Exploration momentum
            % D     --> Critical distance
            %--------------------------------------------------------------
            R0 = obj.Beta(1);
            M = obj.Beta(2);
            D = obj.Beta(3);
        end % getHyperParams
    end % private method signatures
end % classdef