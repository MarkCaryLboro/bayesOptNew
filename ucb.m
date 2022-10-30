classdef ucb < acqFcn
    % Implement the upper confidence bound acquisition function

    properties ( Constant = true )
        FcnName         acqFcnType   = acqFcnType( "ucb" )
    end % Abstract constant properties    
    
    properties ( SetAccess = protected )
        Scale (1,1) double { mustBePositive( Scale ), mustBeReal( Scale ) } = 1
        Beta  (1,1) double { mustBeGreaterThan( Beta, 0 ), ...
                             mustBeLessThan( Beta, 1 ) } = 0.02
    end % protected properties

    methods
        function obj = ucb( ModelObj )
            %--------------------------------------------------------------
            % Upper Confidence Bound (ucb) class constructor
            %
            % obj = ucb( ModelObj );
            %
            % Input Arguments:
            %
            % ModelObj  --> Surrogate model for system
            %--------------------------------------------------------------
            arguments
                ModelObj (1,1)          { mustBeNonempty( ModelObj ) }
            end
            obj.ModelObj = ModelObj;
        end % constructor

        function Fcn = evalFcn( obj, X, Beta )
            %--------------------------------------------------------------
            % Evaluate the UCB acquisition function at the location
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
                obj (1,1)   ucb       {mustBeNonempty( obj )}
                X   (:,:)   double    {mustBeNonempty( X )}
                Beta  (1,1) double    { mustBeGreaterThanOrEqual( Beta, 0),...
                                        mustBeLessThanOrEqual( Beta, 1 )} = obj.Beta
            end
            obj = obj.setBeta( Beta );
            [ Mu, Sigma ] = obj.ModelObj.predict( X );
            Fcn = Mu + sqrt( Beta ) * Sigma;
            Fcn = -Fcn;
        end % evalFcn

        function obj = setBeta( obj, Beta )
            %--------------------------------------------------------------
            % Set the hyperparameter
            %
            % obj = obj.setBeta( Beta )
            %
            % Input Arguments:
            %
            % Beta  --> (double) hyper-parameter value ( 0 <= Beta <= 1 )
            %--------------------------------------------------------------
            arguments
                obj  (1,1)   ucb      {mustBeNonempty( obj )}
                Beta (1,1)   double   { mustBeGreaterThan( Beta, 0),...
                                        mustBeLessThan( Beta, 1 )} = 0.01   
            end
            obj.Beta = Beta;
        end % setBeta

        function obj = setScale( obj, Scale )
            %--------------------------------------------------------------
            % Set the gamma distribution scale parameter
            %
            % obj = obj.setScale( Scale )
            %
            % Input Arguments:
            %
            % Scale --> (double) hyper-parameter value
            %--------------------------------------------------------------
            arguments
                obj  (1,1)   ucb     { mustBeNonempty( obj ) }
                Scale (1,1)  double  { mustBePositive( Scale ),...
                                       mustBeReal( Scale ) } = 1   
            end
            obj.Scale = Scale;
        end % setScale

        function B = sampleGamma( obj, T )
            %--------------------------------------------------------------
            % Sample the gamma distribution to estimate the new Beta value
            %
            % B = obj.sampleGamma( T );
            %
            % Input Arguments:
            %
            % T --> (double) number of iterations
            %--------------------------------------------------------------
            arguments
                obj (1,1)   ucb      { mustBeNonempty( obj ) }
                T   (1,1)   double   { mustBeNonempty( T ),...
                                       mustBePositive( T ) }
            end
            Kt = log( ( T^2 + 1) / sqrt( 2 * pi ) ) / log( 1 + obj.Scale / 2 );
            B = gamrnd( Kt, obj.Scale );
            B = gampdf( B, Kt, obj.Scale );
        end % sampleGamma
    end % Constructor and ordinary methods
end % classdef