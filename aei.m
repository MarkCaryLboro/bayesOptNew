classdef aei < acqFcn
    % Implement the eadaptive xpected improvement acquisition function

    properties ( Constant = true )
        FcnName         acqFcnType   = acqFcnType( "aei" )
    end % Abstract constant properties
    
    properties ( SetAccess = protected )
        Beta                                                                = 1.10
    end

    properties ( SetAccess = protected )
        MnExpRate  (1,:)    double                                          = 0.05
        MxExpRate  (1,:)    double                                          = 10.00
        Alpha      (1,1)    double { mustBeGreaterThanOrEqual( Alpha, 0 ) } = 0.5
        Delta      (1,1)    double { mustBePositive( Delta ) }              = 0.05
        R0         (1,1)    double { mustBePositive( R0 ) }                 = 1.00
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
            % ModelObj --> Surrogate model for system
            % Beta     --> (double) Hyperparameters [ R0, M, Delta, Alpha]
            %--------------------------------------------------------------
            arguments
                ModelObj (1,1)         { mustBeNonempty( ModelObj ) }
                Beta     (1,4) double                                       = [ 1.0, 1.25, 0.05, 0.5 ]      % [ R0, M, D, Alpha ]
            end
            obj.ModelObj = ModelObj;
            %--------------------------------------------------------------
            % Configure the algorithm properties
            %--------------------------------------------------------------
            obj = obj.setInitialExplorationRate( Beta( 1 ) );
            obj = obj.setBeta( Beta( 2 ) );
            obj = obj.setCriticalDistThresh( Beta( 3 ) );
            obj = obj.setPenaltyGain( Beta( 4 ) );
        end % constructor

        function obj = setBeta( obj, M )
            %--------------------------------------------------------------
            % Set the exploration rate multiplier ( M ). The larger M then
            % the more exploratory the search. If M is too small then the
            % algorithm may converge prematurely to a local extremum.
            %
            % obj = obj.setBeta( M );
            %
            % Input Arguments:
            %
            % M --> (double) Exploration rate multiplier (M >= 1)
            %--------------------------------------------------------------
            arguments
                obj     (1,1)  aei    { mustBeNonempty( obj ) }
                M       (1,1)  double { mustBeGreaterThanOrEqual( M, 1 ) } = 1.05
            end
            obj.Beta = M;
        end % setBeta

        function obj = setInitialExplorationRate( obj, R0 )
            %--------------------------------------------------------------
            % Set the initial exploration rate (R0). 
            %
            % Input Arguments:
            %
            % R0 --> (double) Initial exploration rate ( R0 > 0 )
            %--------------------------------------------------------------
            arguments
                obj        (1,1)  aei    { mustBeNonempty( obj ) }
                R0         (1,1)  double { mustBePositive( R0 ) } = 0.05
            end
            obj.R0 = R0;
        end % setInitialExplorationRate

        function obj = setCriticalDistThresh( obj, Delta )
            %--------------------------------------------------------------
            % Set the critical distance threshold which forces an update of
            % the exploration rate (Delta). The larger Delta the more
            % likely the algorithm will become trapped at a local
            % extremum. If too small, convergence may be very slow.
            %
            % Input Arguments:
            %
            % Alpha --> (double) Gain factor for penalty function
            %--------------------------------------------------------------
            arguments
                obj        (1,1)  aei    { mustBeNonempty( obj ) }
                Delta      (1,1)  double { mustBePositive( Delta ) } = 0.001;
            end
            obj.Delta = Delta;
        end % setCriticalDistThresh

        function obj = setPenaltyGain( obj, Alpha )
            %--------------------------------------------------------------
            % Set the penalty gain factor (Alpha). The larger Alpha  then
            % the more exploratory the search. Convergence will be poor if
            % Alpha is too large. Set to zero to disable the penalty.
            %
            % Input Arguments:
            %
            % Alpha --> (double) Gain factor for penalty function
            %--------------------------------------------------------------
            arguments
                obj   (1,1)  aei       { mustBeNonempty( obj ) }
                Alpha (1,1)  double    { mustBeGreaterThanOrEqual( Alpha, 0 ) } = 0.0
            end
            obj.Alpha = Alpha;
        end % setPenaltyGain

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
            % Evaluate the AEI acquisition function at the location
            % specified
            %
            % Fcn = obj.evalFcn( X, Beta );
            %
            % Input Arguments:
            %
            % X     --> Points to evaluate the acquisition function
            % Beta  --> Hyperparameter value
            %--------------------------------------------------------------
            arguments
                obj    (1,1)   aei       {mustBeNonempty( obj )}
                X      (:,:)   double    {mustBeNonempty( X )}
                Beta   (1,:)   double                                       = obj.Beta
            end
            %--------------------------------------------------------------
            % Need a persitent variable for the rate update
            %--------------------------------------------------------------
            persistent Xlast Pen
            if isempty( Xlast ) || obj.Reset
                obj.Reset = false;
                Xlast = X;
                obj.Count = 0;
                Pen = 0;
            end
            obj = obj.setBeta( Beta );
            %--------------------------------------------------------------
            % Threshold test to increase exploration rate
            %--------------------------------------------------------------
            UpdateXlast = ( ( norm( X - Xlast, 1 ) ./ norm( X, 1 ) ) > 0.001 );
            if ( UpdateXlast &&...
                    ( ( norm( X - Xlast, 1 ) / norm( Xlast, 1 ) ) < obj.Delta ) )
                %----------------------------------------------------------
                % Increase the exploration rate
                %----------------------------------------------------------
                obj.Count = obj.Count + 1;
            end
            %--------------------------------------------------------------
            % Calculate the penalty
            %--------------------------------------------------------------
            if UpdateXlast
                Pen = obj.Alpha / norm( X - Xlast, 2 );
            end
            if ( isnan( Pen ) || isinf( Pen ) )
                %----------------------------------------------------------
                % Trap potential error on first pass
                %----------------------------------------------------------
                Pen = 0;
            end
            %--------------------------------------------------------------
            % Calculate the acquisition function
            %--------------------------------------------------------------
            [ Z, Mu, Sigma ] = obj.calcZscore( X );
            Zpdf = normpdf( Z );
            Zcdf = normcdf( Z );
            Explore = obj.ExpMult * obj.R * Sigma .* Zpdf;
            if obj.Problem
                %----------------------------------------------------------
                % Maximisation problem
                %----------------------------------------------------------
                Exploit = ( Mu - obj.Fmax ) .* Zcdf;
                Fcn = Exploit + Explore - Pen;                              % add exploration bonus
                Fcn = -Fcn;                                                 % necessary as fmincon only minimises
            else
                %----------------------------------------------------------
                % Minimisation problem
                %----------------------------------------------------------
                Exploit = ( Mu - obj.Fmin ) .* Zcdf;
                Fcn = Exploit - Explore + Pen;                              % subtract exploration bonus
            end 
            %--------------------------------------------------------------
            % Update Xlast if indicated
            %--------------------------------------------------------------
            if UpdateXlast
                Xlast = X;
            end
        end % evalFcn

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
            if ( Sigma == 0 ) 
                Z = zeros( size( Mu ) );
            elseif matches( string( obj.Problem ), "Maximisation")
                Z = ( Mu - obj.ModelObj.Fmax ) ./ Sigma;
            else
                Z = ( Mu - obj.ModelObj.Fmin ) ./ Sigma;
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
            R = obj.R0 * obj.Beta .^ double( obj.Count );
            %--------------------------------------------------------------
            % Apply clips
            %--------------------------------------------------------------
            R = max( [ R, obj.MnExpRate] );
            R = min( [ R, obj.MxExpRate ] );
        end % get.R
    end % Set/Get methods

    methods ( Access = private )
    end % private method signatures
end % classdef