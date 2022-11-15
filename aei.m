classdef aei < acqFcn
    % Implement the eadaptive xpected improvement acquisition function

    properties ( Constant = true )
        FcnName         acqFcnType   = acqFcnType( "aei" )
    end % Abstract constant properties
    
    properties ( SetAccess = protected )
        Beta  
    end

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
                Beta     (1,4)  double  = [ 0.5, 1, 1.1, 0.25 ]
            end
            obj.Beta = Beta;
            obj.ModelObj = ModelObj;
        end % constructor

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
                Beta    (1,4) double      = obj.Beta
            end
            obj = obj.setBeta( Beta );
            [ A, R0, M, D ] = obj.getHyperParams();
            if ( norm( X - Xstar, 2 ) < D )
                %----------------------------------------------------------
                % Increase the exploration rate
                %----------------------------------------------------------
                R = R0 * M;
                obj = obj.setBeta( [ A, R, M, D ] );
            else
                R = R0;
            end
            P = obj.evalPenalty( X );
            [ Z, Mu, Sigma ] = obj.calcZscore( X );
            Zpdf = normpdf( Z );
            Zcdf = normcdf( Z );
            Fcn = -( Mu - obj.Fmax ) .* Zcdf - R * Sigma .* Zpdf + P;
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
                Beta (1,4)      double  = [ 0.5, 1, 1.1, 0.25 ]
            end
            Ok = Beta >= [ sqrt(eps), 1, 1, sqrt(eps) ];
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
                obj (1,1)   aei      { mustBeNonempty( obj ) }
                X   (:,:)   double  { mustBeNonempty( X ) }                
            end
            [ Mu, Sigma ] = obj.ModelObj.predict( X );
            try
                Z = ( Mu - obj.ModelObj.Fmax ) ./ Sigma;
            catch 
                Z = zeros( size( Mu ) );
            end
        end % calcZscore
    end % ordinary method signatures

    methods ( Access = private )
        function [ A, R, M, D] = getHyperParams( obj )
            %--------------------------------------------------------------
            % Return the hyperparameters
            %
            % [ A, R, M, D] = obj.getHyperParams();
            %
            % Output Arguments:
            %
            % A     --> Penalty function gain
            % R     --> Exploration rate
            % M     --> Exploration momentum
            % D     --> Critical distance
            %--------------------------------------------------------------
            A = obj.Beta(1);
            R = obj.Beta(2);
            M = obj.Beta(3);
            D = obj.Beta(4);
        end % getHyperParams

        function P = evalPenalty( obj, X )
            %--------------------------------------------------------------
            % Evaluate the penalty function: inverse squared distance!
            %
            % P = obj.evalPenalty( X );
            %
            % Input Arguments:
            %
            % X       --> Points to evaluate the acquisition function
            %--------------------------------------------------------------
            A = obj.getHyperParams();
            %--------------------------------------------------------------
            % Carry out the distance comparison in coded units
            %--------------------------------------------------------------
%             Xc = obj.ModelObj.code( X );
%             Xstar = obj.ModelObj.code( obj.BestX );
            %--------------------------------------------------------------
            % Calculate Euclidean distance to best known result
            %--------------------------------------------------------------
            P = vecnorm( X - obj.BestX, 2, 2 );
            if P ~= 0
                %----------------------------------------------------------
                % Normal case
                %----------------------------------------------------------
                P = A ./ P.^2;
            else
                %----------------------------------------------------------
                % Trap infinity penalty case
                %----------------------------------------------------------
                P = A ./ 0.001;
            end
        end % evalPenalty
    end % private method signatures
end % classdef