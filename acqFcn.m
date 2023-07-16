classdef ( Abstract = true) acqFcn < handle
    % Interface to Bayesian Optimisation acquisition functions

    properties ( Constant = true, Abstract = true)
        FcnName         acqFcnType
    end % Abstract constant properties

    properties ( SetAccess = protected, Abstract = true )
        Beta                                                                % Hyperparameter
    end

    properties ( SetAccess = protected )
        ModelObj (1,1)                                                      % Surrogate model object
        BestX    (1,:)  double  = -inf
        Problem  (1,1)  logical = true
    end % protected properties

    properties ( SetAccess = protected, Dependent )
        Data     double                                                     % Current input training data
        Response double                                                     % Current training response data
    end % dependent properties

    properties ( Access = protected, Dependent )
        Fmax     double                                                     % Maximum function evaluation
        Nvar     double                                                     % Number of input variables               
    end % protected dependent properties

    methods ( Abstract = true )
        Y = evalFcn( obj, X )                                               % Evaluate the acquisition function at the points supplied
        obj = setBeta( obj, Beta )                                          % Set the hyper-parameter value
    end % Abstract method signatures

    methods
        function obj = setModelType( obj, M )
            %--------------------------------------------------------------
            % Set the problem type property ( Problem ).
            %
            % obj = obj.setModelType( M );
            %
            % Input Arguments:
            %
            % M --> (logical) set to true for a maximisation problem, and
            %                 false for a minimisation problem.
            %--------------------------------------------------------------
            arguments
                obj (1,:)   
                M   (1,1) logical   = true
            end
            obj.Problem = M;
        end % setModelType

        function obj = addFcnSample2Data( obj, Xnew, Ynew )
            %--------------------------------------------------------------
            % Add a new function sample to the training data for the
            % surrogte model and retrain it. 
            %
            % obj = obj.addFcnSample( Xnew, Ynew );
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
            obj.ModelObj = obj.ModelObj.updateModel( Xnew, Ynew );
        end % addFcnSample2Data
    end % ordinary method signatures

    methods
        function D = get.Data( obj )
            % Return the current training data inputs
            D = obj.ModelObj.X;
        end % get.Data

        function Y = get.Response( obj )
            % Return the current training data response vector
            Y = obj.ModelObj.Y;
        end % get.Response

        function N = get.Nvar( obj )
            % Return number of variables
            N = size( obj.Data, 2 );
        end % get.Nvar

        function Fmax = get.Fmax( obj )
            % Return best value of acquisition function
            try
                Fmax = obj.ModelObj.Fmax;
            catch
                Fmax = NaN;
            end
        end % get.Fmax
    end % get/set method signatures
end % classdef