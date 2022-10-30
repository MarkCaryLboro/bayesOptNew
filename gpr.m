classdef gpr < surrogateModel
    % gaussian regression process class

    properties ( Constant = true )
        ModelType   string    = "GP"
    end % Abstract constnat properties

    properties ( SetAccess = protected )
        ModelObj    RegressionGP
        Kernel      kernels         = kernels( "ARDsquaredExponential" )
        PredMethod  gprPredMethod   = gprPredMethod( "exact" )
        FitMethod   gprFitMethod    = gprFitMethod( "exact" )
    end % protected and abstract property definitions

    properties ( SetAccess = protected, Dependent )
        LenScale    double                                                  % Length scales vector
        SigmaF      double                                                  % Noise covariance function
    end % dependent properties

    properties ( Access = protected, Dependent )
        Cov         double                                                  % Kernel matrix for the training data
    end % protected dependent properties

    methods
        function obj = gpr( X, Y )
            %--------------------------------------------------------------
            % Class constructor
            %
            % obj = gpr( X, Y );
            %
            % Input Arguments:
            %
            % X --> (double) Input data matrix
            % Y --> (double) Response data matrix
            %--------------------------------------------------------------
            arguments
                X   double     = [];
                Y   double     = [];
            end
            obj = obj.setTrainingData( X, Y );
        end % gpr
    end % constructor method signature    

    methods
        function obj = setPredMethod( obj, Method )
            %--------------------------------------------------------------
            % Set the method used to make predictions
            %
            % obj = obj.setPredMethod( Method );
            %
            % Input Arguments:
            %
            % Name  --> (string) Supported prediction method:
            %           exact {default}
            %           sd           
            %           sr              
            %           fic             
            %           bcd 
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   gpr
                Method  (1,1)   string  = "exact"
            end

            try
                M = gprPredMethod( Method );
                obj.PredMethod = M;
            catch me
                error( me.identifier , me.message );
            end        
        end % setPredMethod

        function obj = setFitMethod( obj, Method )
            %--------------------------------------------------------------
            % Set the method used to fit the data
            %
            % obj = obj.setFitMethod( Method );
            %
            % Input Arguments:
            %
            % Name  --> (string) Supported prediction method:
            %           exact {default}
            %           sd           
            %           sr              
            %           fic             
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   gpr
                Method  (1,1)   string  = "exact"
            end

            try
                M = gprFitMethod( Method );
                obj.FitMethod = M;
            catch me
                error( me.identifier , me.message );
            end        
        end % setFitMethod

        function obj = setKernel( obj, Name )
            %--------------------------------------------------------------
            % Set the gaussian covariance function kernel
            %
            % obj = obj.setKernel( Name );
            %
            % Input Arguments:
            %
            % Name  --> (string) Supported kernel function:
            %           ARDsquaredExponential {default}
            %           ARDexponential           
            %           ARDmatern32              
            %           ARDmatern52             
            %--------------------------------------------------------------
            arguments
                obj     (1,1)   gpr
                Name    (1,1)   string  = "ARDsquaredExponential"
            end
            try
                K = kernels( Name );
                obj.Kernel = K;
            catch me
                error( me.identifier , me.message );
            end
        end % setKernel

        function obj = trainModel( obj, varargin )
            %--------------------------------------------------------------
            % Train the gpr model.
            %
            % obj = trainModel( obj, Name1, Value1, ..., Name#, Value#);
            %
            % This method makes use of the RegressionGP class. A
            % RegressionGP object is composited in the ModelObj property.
            %
            % The kernel function, fit and prediction methods can be set 
            % using the setKernel and setPredMethod methods. However, all 
            % other options can be set through the ( Name, Value ) pair 
            % list for fitrgp().
            %--------------------------------------------------------------
            assert( obj.DataOk,...
                    "X-predictor name vector must have %3.0f entries",...
                    obj.N );            
            obj.ModelObj = fitrgp( obj.Xc, obj.Y,...
                    "FitMethod", string( obj.FitMethod ),...
                    "PredictMethod", string( obj.PredMethod ),...
                    "KernelFunction", string( obj.Kernel ),...
                    varargin{:});
            obj.Trained = true;
        end % trainModel

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
                obj   (1,1)     gpr
                Xnew            double           = obj.X     
                Alpha           double  = 0.05
            end
            assert( obj.Trained, 'Must first train the model using the "trainModel" method before predictions can be made')
            Xnew = obj.code( Xnew );
            [ Ypred, Ysd, Yint ] = predict( obj.ModelObj, Xnew, 'Alpha',...
                                            Alpha );
        end % predict

        function K = genK( obj, X, Xref )
            %--------------------------------------------------------------
            % Generate the kernel matrix
            %
            % K = obj.genK( X );
            %
            % Input Arguments:
            %
            % X --> Input data matrix.
            %--------------------------------------------------------------
            arguments 
                obj     (1,1)       gpr
                X       (:,:)       double       = obj.X
                Xref    (:,:)       double       = obj.X
            end
            Xc = obj.code( X );
            Xc_ref = obj.code( Xref );
            switch obj.Kernel
                case "ARDsquaredExponential"
                    K = obj.ardSquaredExponential( Xc, Xc_ref );
                case "ARDexponential"
                    K = obj.ardExponential( Xc, Xc_ref );
                case "ARDmatern32"
                    K = ardMatern32( Xc, Xc_ref );
                case "ARDmatern52"
                    K = ardMatern52( Xc, Xc_ref );
                otherwise
                    error( 'Unknown Kernel Type "%s', obj.Kernel );
            end
        end % genK

        function S = sigma( obj, X )
            %--------------------------------------------------------------
            % Covariance matrix for predictions
            %
            % S = obj.sigma( X );
            %
            % Input Arguments:
            %
            % X --> Location in the Gaussian field to calculate sigma.
            %--------------------------------------------------------------
            arguments 
                obj (1,1)       gpr
                X   (:,:)       double       = obj.X
            end
            %--------------------------------------------------------------
            % Generate required kernel matrices
            %--------------------------------------------------------------
            K = obj.genK( obj.X );
            Kstar = obj.genK( X, X );
            KsK = obj.genK( X, obj.X );
            %--------------------------------------------------------------
            % Calculate the corresponding covarince matrix
            %--------------------------------------------------------------
            S = ( K + obj.SigmaF^2 * eye( obj.NumPoints ) ) \ ...
                      eye( obj.NumPoints );
            S = KsK.' * S * KsK;
            S = Kstar -  S;
        end % sigma
    end % ordinary methods signatures

    methods
        function K = get.Cov( obj )
            % Return the covariance matrix for the training data
            K = obj.sigma( obj.X );
        end

        function L = get.LenScale( obj )
            % Return the lengthscale parameters
            L = obj.ModelObj.KernelInformation.KernelParameters;
            L = L( 1:( end - 1 ) );
        end

        function SF = get.SigmaF( obj )
            % Return the estimated noise scale
            SF = obj.ModelObj.KernelInformation.KernelParameters( end );
        end
    end % get/set methods

    methods (  Access = private )
        function K = ardSquaredExponential( obj, Dc, Xc )
            %--------------------------------------------------------------
            % Calculate the Kernel function, for the data supplied
            %
            % K = obj.ardSquaredExponential( Dc );
            %
            % Input Arguments:
            %
            % Dc    -> (Nxp) matrix of coded input data
            % Xc    -> (Mxp) matrix of coded reference data
            %--------------------------------------------------------------
            R = obj.genScaledSqrdResMatrix( Dc, Xc );
            K = obj.SigmaF^2 .* exp( -0.5 * R );
        end % ardSquaredExponential

        function K = ardExponential( obj, Dc, Xc )
            %--------------------------------------------------------------
            % Calculate the Kernel function, for the data supplied
            %
            % K = obj.ardExponential( Dc );
            %
            % Input Arguments:
            %
            % Dc    -> (Nxp) matrix of coded input data
            % Xc    -> (Mxp) matrix of coded reference data
            %--------------------------------------------------------------
            R = obj.genScaledSqrdResMatrix( Dc, Xc );
            K = obj.SigmaF^2 .* exp( -sqrt( R ) );
        end % ardExponential

        function K = ardMatern32( obj, Dc, Xc )
            %--------------------------------------------------------------
            % Calculate the Kernel function, for the data supplied
            %
            % K = obj.rdMatern32( Dc );
            %
            % Input Arguments:
            %
            % Dc    -> (Nxp) matrix of coded input data
            % Xc    -> (Mxp) matrix of coded reference data
            %--------------------------------------------------------------   
            R = sqrt( obj.genScaledSqrdResMatrix( Dc, Xc ) );
            K = obj.SigmaF^2 * ( 1 + sqrt( 3 ) * R ) .* exp( -sqrt( 3 ) * R );            
        end % ardMatern32

        function K = ardMatern52( obj, Dc, Xc )
            %--------------------------------------------------------------
            % Calculate the Kernel function, for the data supplied
            %
            % K = obj.ardMatern52( Dc );
            %
            % Input Arguments:
            %
            % Dc    -> (Nxp) matrix of coded input data
            % Xc    -> (Mxp) matrix of coded reference data
            %--------------------------------------------------------------   
            R = sqrt( obj.genScaledSqrdResMatrix( Dc, Xc ) );
            K = obj.SigmaF^2 * ( 1 + sqrt( 5 ) * R + ( 5 / 3 ) * R.^2 )...
                .* exp( -sqrt( 5 ) * R );            
        end % ardMatern52

        function R = genScaledSqrdResMatrix( obj, Dc, Xc )
            %--------------------------------------------------------------
            % Generate sum sqaured residual matrix for data supplied
            %
            % R = obj.genScaledSqrdResMatrix( Dc );
            %
            % Input Arguments:
            %
            % Dc    -> (Nxp) matrix of coded input data
            % Xc    -> (Mxp) matrix of coded reference data
            %--------------------------------------------------------------   
            D = diag( 1 ./ obj.LenScale );
            Dc = Dc * D;
            Pxc = size( Xc, 1 );
            Pdc = size( Dc, 1 );
            R = zeros( Pxc, Pdc );
            for Q = 1:Pdc
                D_q = Dc( Q, : );
                R( :, Q ) = sum( ( Xc - D_q ).^2, 2 );
            end
        end % genScaledSqrdResMatrix
    end % Private methods signatures

    methods ( Static = true, Access = private )
    end % Static and private methods
end % classdef