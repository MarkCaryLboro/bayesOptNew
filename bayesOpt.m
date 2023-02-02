classdef bayesOpt < handle
    % Bayesian optimisation class

    events
        UPDATE                                                              % Run newly requested simulation                                                   
    end

    properties ( SetAccess = protected, Dependent  )
        AcqFcn   (1,1)   string                                             % Acquisition function type
        Model    (1,1)   string                                             % Surrogate model type
        X        (:,:)   double                                             % Current sample input data
        Y        (:,1)   double                                             % Current function query data
        Xnext    (1,:)   double                                             % Next point to query
        ModelObj (1,1)                                                      % Surrogate model object
        HyperPar (1,1)   double                                             % Hyper-parameter
        Xbest    (1,:)   double                                             % Best known input configuration
        Ybest    (1,1)   double                                             % Best known function query
        Xlo      (1,:)   double                                             % Lower bound for x-data (for data coding)
        Xhi      (1,:)   double                                             % Upper bound for x-data (for data coding)
    end % dependent properties

    properties ( SetAccess = protected )
        AcqObj   (1,1)                                                      % Acquisition function object
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
            %               either {"ucb"} or "ei".
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
            obj.AcqObj.setBestX( M.Xmax );
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
                A   (1,:)   double      { mustBeNonempty( A ) }
                B   (1,:)   double      { mustBeNonempty( B ) }
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
                Args.lb        (1,:) double                
                Args.ub        (1,:) double                
                Args.options   (1,1) optim.options.Fmincon = optimoptions( "fmincon")
                Args.nonlincon (1,1) function_handle       
                Args.Aineq     (:,:) double                
                Args.bineq     (:,1) double                
                Args.Aeq       (:,:) double                
                Args.beq       (:,1) double                
            end
            %--------------------------------------------------------------
            % Parse the optional arguments
            %--------------------------------------------------------------
            Names = [ "lb", "ub", "nonlincon", "Aineq", "binq", "Aeq",...
                       "beq", "options" ];
            for Q = 1:numel( Names )
               try
                   Problem.( Names( Q ) ) = Args.( Names( Q ) );
               catch
                   Problem.( Names( Q ) ) = [];
               end
            end
            %--------------------------------------------------------------
            % Set up the optimisation problem
            %--------------------------------------------------------------
            Problem.options.Display = "iter";
            Problem.solver = "fmincon"; 
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
            Problem.objective = @(X)obj.AcqObj.evalFcn( X, B );
            if isinf( obj.Xbest )
                %----------------------------------------------------------
                % Select a starting point
                %----------------------------------------------------------
                [ ~, Idx ] = max( obj.Y );
                Xmax = obj.X( Idx, : );
                obj.AcqObj = obj.AcqObj.setBestX( Xmax );
            end
            Problem.x0 = obj.Xbest;
            BestX = fmincon( Problem );
            obj.AcqObj.setBestX( BestX );
        end % acqFcnMaxTemplate

        function broadcastUpdate( obj )
            %--------------------------------------------------------------
            % Broadcast the UPDATE message to generate a new function query
            %
            % obj.broadcastUpdate();
            %--------------------------------------------------------------
            notify( obj, "UPDATE");
        end % 
    end % Concrete ordinary methods

    methods
        function Xnext = get.Xnext( obj )
            Xnext = obj.AcqObj.BestX;
        end % get.Xnext
        
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
            Y = max( obj.Y );
        end % get.Ybest

        function X = get.Xbest( obj )
            % Return the input configuration from the sample pool
            % corresponding to the best function query
            [ ~, Idx ] = max( obj.Y );
            X = obj.X( Idx, : );
        end % get.Xbest

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
end % classdef