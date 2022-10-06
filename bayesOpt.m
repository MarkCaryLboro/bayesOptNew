classdef bayesOpt 
    % Bayesian optimisation class

    properties ( SetAccess = protected, Dependent  )
        AcqFcn   (1,1)   string                                             % Acquisition function type
        Model    (1,1)   string                                             % Surrogate model type
        X        (:,:)   double                                             % Current sample input data
        Y        (:,1)   double                                             % Current function query data
        Xnext    (1,:)   double                                             % Next point to query
        ModelObj (1,1)                                                      % Surrogate model object
    end % dependent properties

    properties ( SetAccess = protected )
        AcqObj   (1,1)                                                      % Acquisition function object
    end % Portected properties
    
    methods
        function obj = bayesOpt( AcqFcnObj )
            %--------------------------------------------------------------
            % Bayesian optimisation class constructor
            %
            % obj = bayesOpt( AcqFcnObj );
            %
            % Input Arguments:
            %
            % AcqFcnObj     --> Acquisition function object, either an "ei"
            %                   or "ucb" object.
            %--------------------------------------------------------------
            arguments
                AcqFcnObj (1,1)    { mustBeAcqFcn( AcqFcnObj ) }
            end
            obj.AcqObj = AcqFcnObj;
        end % constructor

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
            %==============================================================
            % To Do: Generalise the acquisition function object to have a
            % hyper-parameter Dependent property that is an abstract member
            %==============================================================
            switch lower( class( obj.AcqObj ) )
                case 'ei'
                    Problem.objective = @(X)obj.AcqObj.evalFcn( X,...
                                                       obj.AcqObj.Xi );
                case 'ucb'
                    B = obj.AcqObj.sampleGamma( numel (obj.Y) );
                    Problem.objective = @(X)obj.AcqObj.evalFcn( X, B );
                otherwise
            end
            if isinf( obj.AcqObj.BestX )
                %----------------------------------------------------------
                % Select a strating point
                %----------------------------------------------------------
                [ ~, Idx ] = max( obj.AcqObj.ModelObj.predict( obj.X ) );
                Xmax = obj.AcqObj.Data( Idx, : );
                obj.AcqObj = obj.AcqObj.setBestX( Xmax );
            end
            Problem.x0 = obj.AcqObj.BestX;
            BestX = fmincon( Problem );
            obj.AcqObj = obj.AcqObj.setBestX( BestX );
        end % acqFcnMaxTemplate
    end % Concrete ordinary methods

    methods
        function Xnext = get.Xnext( obj )
            Xnext = obj.AcqObj.BestX;
        end % get.Xnext

        function X = get.X( obj )
            % Return the current function query input locations
            X = obj.AcqObj.Data;
        end % get.X

        function Y = get.Y( obj )
            % Return the current function query input locations
            Y = obj.AcqObj.Response;
        end % get.Y

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

function Ok = mustBeAcqFcn( This )
    try
        FcnName = string( This.FcnName );
        Ok = any( contains( FcnName, [ "ei", "ucb" ] ) );
    catch
        Ok = false;
    end
end