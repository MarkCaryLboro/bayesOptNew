classdef bayesOpt < handle
    % Bayesian optimisation class

    properties ( SetAccess = protected, Dependent  )
        AcqFcn   (1,1)   string                                             % Acquisition function type
        Model    (1,1)   string                                             % Surrogate model type
        X        (:,:)   double                                             % Current sample input data
        Y        (:,1)   double                                             % Current function query data
        ModelObj (1,1)                                                      % Surrogate model object
        HyperPar (1,1)   double                                             % Hyper-parameter
        Xbest    (1,:)   double                                             % Best known input configuration
        Ybest    (1,1)   double                                             % Best known function query
        Xlo      (1,:)   double                                             % Lower bound for x-data (for data coding)
        Xhi      (1,:)   double                                             % Upper bound for x-data (for data coding)
        Bidx     (1,1)   double                                             % Pointer to best result
        Problem  (1,1)   optimisationType                                   % Either "maximisation" or "minimisation"
    end % dependent properties

    properties ( SetAccess = protected )
        AcqObj   (1,1)                                                      % Acquisition function object
        Xnext    (1,:)   double                                             % Next point to query
    end % Portected properties
    
    methods
        function obj = bayesOpt( Model, AcqFcn )
            %--------------------------------------------------------------
            % Bayesian optimisation class constructor
            %
            % obj = bayesOpt( Modeltype, AcqFcn );
            %
            % Input Arguments:
            %
            % Model     --> (string) surrogate model type. Must be either
            %               {"gpr"} or "rf".
            % AcqFcn    --> (string) Acquisition function name. Must be 
            %               either {"ucb"}, "aei" or "ei".
            %--------------------------------------------------------------
            arguments
                Model  (1,1)    string    = "gpr"
                AcqFcn (1,1)    string    = "ucb"
            end
            Model = surrogateModelType( Model );
            ModelObj = eval( lower( string( Model ) ) );                    %#ok<NASGU> 
            AcqFcn = acqFcnType( AcqFcn );
            Str = lower( string ( AcqFcn ) );
            Str = strjoin( [ Str, "( ModelObj )"], "" );
            obj.AcqObj = eval( Str );
        end % constructor

        function obj = setHyperPar( obj, Beta )
            %--------------------------------------------------------------
            % Set the hyperparameter value for the acquisition function
            % optimisation problem
            %
            % obj = obj.setHyperPar( Beta );
            %
            % Input Arguments:
            %
            % Beta  --> (double) Hyperparameter value
            %--------------------------------------------------------------
            arguments
                obj     (1,1)           { mustBeNonempty( obj ) }
                Beta    (1,1)   double
            end
            A = obj.AcqObj;
            A.setBeta( Beta );
        end % setHyperPar

        function obj = setProblemTypeState( obj, M )
            %--------------------------------------------------------------
            % Set the problem type property ( MaxFnc ).
            %
            % obj = obj.setProblemTypeState( M );
            %
            % Input Arguments:
            %
            % M --> (string) set to "Maximisation" or "Minimisation" as 
            %                appropriate.
            %--------------------------------------------------------------
            arguments
                obj (1,:) bayesOpt  
                M   (1,1) string          { mustBeMember( M, ...
                                            ["Maximisation",...
                                             "Minimisation"] ) }            = "Maximisation"
            end
            obj.AcqObj.setProblemType( M );
        end % setProblemTypeState

        function obj = setTrainingData( obj, X, Y )
            %--------------------------------------------------------------
            % Set the training data
            %
            % obj = obj.setTrainingData( X, Y );
            %
            % Input Arguments:
            %
            % X --> (double) NxC matrix of inputs
            % Y --> (double) Nx1 matrix of response data
            %--------------------------------------------------------------
            arguments
                obj (1,1)
                X           double
                Y   (:,1)   double
            end
            M = obj.AcqObj.ModelObj;
            %--------------------------------------------------------------
            % This works because the surrogate model class inherits from
            % the handle class...
            %--------------------------------------------------------------
            M.setTrainingData( X, Y );
            M.trainModel;
        end % setTrainingData

        function [ Ypred, Ysd, Yint ] = predict( obj, Xnew, Alpha )
            %--------------------------------------------------------------
            % Model predictions
            %
            % Y = obj.predict( Xnew, Alpha );
            %
            % Input Argumentss:
            %
            % Xnew  --> (double) input data
            % Alpha --> (double) 100(1 - Alpha)% prediction interval
            %
            % Output Arguments:
            %
            % Ypred --> predicted responses
            % Ysd   --> standard deviations
            % Yint  --> prediction interval
            %--------------------------------------------------------------
            arguments
                obj   (1,1)     
                Xnew            double          = obj.AcqObj.ModelObj.X     
                Alpha           double          = 0.05
            end
            [ Ypred, Ysd, Yint ] = obj.AcqObj.ModelObj.predict( Xnew,...
                                            Alpha );
        end % predict

        function obj = addNewQuery( obj, Xnew, Ynew )
            %--------------------------------------------------------------
            % Add a new function query & update surrogate model.
            %
            % obj = obj.addNewQuery( Xnew, Ynew );
            %
            %
            % Input Arguments:
            %
            % Xnew  --> (double) New input data
            % Ynew  --> (double) Function query f(Xnew)
            %--------------------------------------------------------------
            arguments
                obj     (1,1)           { mustBeNonempty( obj ) }
                Xnew    (:,:)   double  { mustBeNonempty( Xnew ) }
                Ynew    (:,1)   double  { mustBeNonempty( Ynew ) }
            end
            obj.AcqObj.addFcnSample2Data( Xnew, Ynew );
        end % addNewQuery

        function obj = conDataCoding( obj, A, B )
            %--------------------------------------------------------------
            % Configure the data coding
            %
            % obj = obj.conDataCoding( A, B );
            %
            % Input Arguments:
            %
            % A     --> (double) Lower bound for x-data 
            % B     --> (double) Upper bound for x-data 
            %--------------------------------------------------------------
            arguments
                obj (1,1)   bayesOpt    { mustBeNonempty( obj ) }
                A   (1,:)   double      { mustBeNonempty( A ) } = min( obj.X )
                B   (1,:)   double      { mustBeNonempty( B ) } = max( obj.X )
            end
            M = obj.ModelObj;
            M.conDataCoding( A, B );
        end % conDataCoding
        
        function obj = acqFcnMaxTemplate( obj, Args )
            %--------------------------------------------------------------
            % Optimise the acquisition function given the current training
            % data. The output is the next place to sample the function.
            %
            % Xmax = obj.acqFcnTemplate( "Name1", Value1, ... );
            %
            % Input Arguments:
            %
            % "Name"    --> Fmincon input argument name
            % Value     --> Corresponding argument value
            %--------------------------------------------------------------
            arguments
                obj            (1,1)            { mustBeNonempty( obj ) }
                Args.lb        (1,:) double                = obj.Xlo          
                Args.ub        (1,:) double                = obj.Xhi
                Args.options   (1,1) optim.options.Fmincon = optimoptions( "fmincon")
                Args.nonlcon   (1,1) function_handle       
                Args.Aineq     (:,:) double                
                Args.bineq     (:,1) double                
                Args.Aeq       (:,:) double                
                Args.beq       (:,1) double                
            end
            %--------------------------------------------------------------
            % Parse the optional arguments
            %--------------------------------------------------------------
            Names = [ "lb", "ub", "nonlcon", "Aineq", "binq", "Aeq",...
                       "beq", "options" ];
            for Q = 1:numel( Names )
               try
                   PROBLEM.( Names( Q ) ) = Args.( Names( Q ) );
               catch
                   PROBLEM.( Names( Q ) ) = [];
               end
            end
            %--------------------------------------------------------------
            % Code the bound constraints
            %--------------------------------------------------------------
            PROBLEM.lb = obj.code( PROBLEM.lb );
            PROBLEM.ub = obj.code( PROBLEM.ub );
            %--------------------------------------------------------------
            % Set up the optimisation problem
            %--------------------------------------------------------------
            PROBLEM.options.Display = "iter";
            PROBLEM.solver = "fmincon"; 
            %--------------------------------------------------------------
            % Set the hyperparameter
            %--------------------------------------------------------------
            switch lower( class( obj.AcqObj ) )
                case 'ucb'
                    B = obj.AcqObj.sampleGamma( numel(obj.Y) );
                otherwise
                    B = obj.HyperPar;
            end
            switch lower( class( obj.AcqObj ) )
                case 'aei'
                    % reset the exploration rate index counter every time
                    % the optimisation is called
                    obj.AcqObj.resetCounter();
            end
            PROBLEM.objective = @(X)obj.AcqObj.evalFcn( X, true, B );
            if isempty( obj.Xnext )
                obj = obj.setXbestAsXnext();
            end
            PROBLEM.x0 = obj.code( obj.Xnext );
            Xnew = fmincon( PROBLEM );
            obj.Xnext = obj.decode( Xnew );
        end % acqFcnMaxTemplate

        function obj = setXbestAsXnext( obj )
            %--------------------------------------------------------------
            % Set the next query to the current best value
            %
            % obj = obj.setXbestAsXnext();
            %--------------------------------------------------------------
            arguments
                obj (1,1) bayesOpt { mustBeNonempty( obj ) }
            end
            obj.Xnext = obj.Xbest;
        end % setXbestAsXnext
    end % Concrete ordinary methods

    methods        
        function B = get.HyperPar( obj )
            B = obj.AcqObj.Beta;
        end % get.Beta
            
        function X = get.X( obj )
            % Return the current function query input locations
            X = obj.AcqObj.Data;
        end % get.X

        function A = get.Xlo( obj )
            % Return the lower data bound for coding
            A = obj.ModelObj.Xlo;
        end % get.Xlo

        function B = get.Xhi( obj )
            % Return the lower data bound for coding
            B = obj.ModelObj.Xhi;
        end % get.Xhi

        function Y = get.Y( obj )
            % Return the current function query input locations
            Y = obj.AcqObj.Response;
        end % get.Y

        function Y = get.Ybest( obj )
            % Return the best query from the sample pool
            switch obj.Problem
                case "Maximisation"
                    Y = max( obj.Y );
                otherwise
                    Y = min( obj.Y );
            end
        end % get.Ybest

        function X = get.Xbest( obj )
            % Return the input configuration from the sample pool
            % corresponding to the best function query
            X = obj.X( obj.Bidx, : );
        end % get.Xbest

        function P = get.Problem( obj )
            % return the problem type
            P = obj.AcqObj.Problem;
        end % get.Problem

        function B = get.Bidx( obj )
            % Return the position of the best result in the data pool
            switch obj.Problem
                case "Maximisation"
                    [ ~, B ] = max( obj.Y );
                otherwise
                    [ ~, B ] = min( obj.Y );
            end
        end % get.Bidx

        function M = get.Model( obj )
            % Return the surrogate model type
            M = obj.AcqObj.ModelObj.ModelType;
        end % get.Model

        function M = get.ModelObj( obj )
            % Return the surrogate model object
            M = obj.AcqObj.ModelObj;
        end % get.ModelObj

        function A = get.AcqFcn( obj )
            % Return the acquisition function nane
            A = obj.AcqObj.FcnName;
        end % get.AcqFcn
    end % get/set methods

    methods ( Access = protected )
        function Xc = code( obj, X )
            %--------------------------------------------------------------
            % Map the input data [a, b] onto the interval [0, 1]
            %
            % Input Arguments:
            %
            % X --> (double) data matrix
            %--------------------------------------------------------------
            M = obj.AcqObj.ModelObj;
            Xc = M.code( X );
        end % code        end % code

        function X = decode( obj, Xc )
            %--------------------------------------------------------------
            % Map the coded input data [-1, 1] onto the interval [a, b]
            %
            % D = obj.decode( Dc );
            %
            % Input Arguments:
            %
            % Dc --> Coded input data
            %--------------------------------------------------------------
            M = obj.AcqObj.ModelObj;
            X = M.decode( Xc );
        end % decode
    end % protected methods
end % classdef