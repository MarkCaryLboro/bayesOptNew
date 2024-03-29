classdef ( Abstract = true ) surrogateModel < handle
    
    properties ( Constant = true, Abstract = true )
        ModelType   string
    end % Abstract constnat properties

    properties ( SetAccess = protected )
        X           double                                                  % Input data
        Y           double                                                  % Response data
        Yname       string                                                  % Response name 
        Xname       string                                                  % Array of predictor names
        Trained     logical         = false                                 % Model trained state flag
        Xunits      string                                                  % Units for input variables
        Yunits      string                                                  % Units for response variable
        Xlo         double                                                  % Lower bound for x-data (for data coding)
        Xhi         double                                                  % Upper bound for x-data (for data coding)
    end % protected and abstract properties

    properties ( Access = protected, Dependent, Abstract = true )
        Cov         double                                                  % Kernel matrix for the training data
    end

    properties ( Access = protected, Dependent )
        Xc          double                                                  % Coded input data [ a, b ] --> [ -1, 1 ]
    end % protected & dependent properties

    properties ( SetAccess = protected, Dependent )
        DataOk      logical                                                 % True if data dimensions are consistent
        N           int8                                                    % Number of variables
        NumPoints   int64                                                   % Number of data points
        Fmax        double                                                  % Maximum function query discovered so far
        Xmax        double                                                  % Location of maximum function query
        Fmin        double                                                  % Minimum function query discovered so far
        Xmin        double                                                  % Location of minimum function query
    end % dependent properties

    methods ( Abstract = true )
        obj = trainModel( obj, varargin )
        Y = predict( obj, X, varargin )
        S = sigma( obj, X, varargin )
    end 

    methods
        function obj = updateModel( obj, Xnew, Ynew )
            %--------------------------------------------------------------
            % Update the model based on additional training data
            %
            % obj = obj.updateModel( Xnew, Ynew );
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
            Xnew = [ obj.X ; Xnew ];
            Ynew = [ obj.Y; Ynew];
            obj = obj.setTrainingData( Xnew, Ynew );
            obj = obj.trainModel();
        end % updateModel

        function obj = setTrainingData( obj, X, Y )
            %--------------------------------------------------------------
            % Set the training data
            %
            % obj = obj.setTrainingData( X, Y );
            %
            % Input Arguments:
            %
            % X --> (double) NxC matrix of input data
            % Y --> (double) Nx1 matrix of response data
            %--------------------------------------------------------------
            arguments
                obj (1,1)
                X   (:,:)   double
                Y   (:,1)   double
            end
            obj.Trained = false;
            obj.X = X;
            obj.Y = Y;
        end % setTrainingData

        function D = decode( obj, Dc )
            %--------------------------------------------------------------
            % Map the coded input data [-1, 1] onto the interval [a, b]
            %
            % D = obj.decode( Dc );
            %
            % Input Arguments:
            %
            % Dc --> Coded input data
            %--------------------------------------------------------------
            [ A, B ] = obj.dataCodingInfo( );
            D = Dc .* ( B - A ) + A;
        end % decode

        function obj = setVarUnits( obj, Xunits, Yunits )
            %--------------------------------------------------------------
            % Set the input and response unit names.
            %
            % obj = obj.setVarUnits( Xunits, Yunits );
            %
            % Input Arguments:
            %
            % Xunits    --> (1xp) (string) vector of x-unit names
            % Yunits    --> (1x1) (string) response variable units
            %--------------------------------------------------------------
            arguments
                obj     (1,1)             
                Xunits  (1,:)   string
                Yunits  (1,1)   string
            end
            obj.Xunits = Xunits;
            obj.Yunits = Yunits;
        end % setVarUnits

        function obj = setVarNames( obj, Xnames, Yname )
            %--------------------------------------------------------------
            % Set the labels for the predictor and response variables
            %
            % obj = obj.setVarNames( Xnames, Yname );
            %
            % Input Arguments:
            %
            % Xnames    --> (string) list of predictor variable names
            % Yname     --> (string) response variable name
            %--------------------------------------------------------------
            arguments
                obj     (1,1)
                Xnames  (1,:)   string  { mustBeNonempty(Xnames) }
                Yname   (1,1)   string    = "Y"
            end
            Ok = ( numel( Xnames ) == obj.N );
            assert( Ok, "X-predictor name vector must have %3.0f entries",...
                    obj.N );
            obj.Xname = Xnames;
            obj = obj.setYname( Yname );
        end % setVarNames

        function obj = setYname( obj, Yname )
            %--------------------------------------------------------------
            % Set the response variable name
            %
            % obj = obj.setYname( Yname );
            %
            % Input Arguments:
            %
            % Yname     --> (string) response variable name
            %--------------------------------------------------------------
            arguments
                obj     (1,1)
                Yname   (1,1)   string    = "Y"
            end
            obj.Yname = Yname;
        end % setYname

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
                obj     (1,1)               { mustBeNonempty( obj ) }
                A       (1,:)   double
                B       (1,:)   double
            end
            try
                obj = obj.setXlo( A );
            catch
                % apply default
                obj = obj.setXlo();
            end
            try
                obj = obj.setXhi( B );
            catch
                % apply default
                obj = obj.setXhi();
            end
        end % conDataCoding

        function Xc = code( obj, X )
            %--------------------------------------------------------------
            % Code data onto the interval [0,1]
            %
            % Input Arguments:
            %
            % X --> (double) data matrix
            %--------------------------------------------------------------
            [ A, B ] = obj.dataCodingInfo();
            Xc = ( X - A ) ./ ( B - A );
        end % code
    end % ordinary method signatures

    methods ( Access = protected )
        function [ A, B, M ] = dataCodingInfo( obj )
            %--------------------------------------------------------------
            % Return data coding information
            %
            % [ A, B, M ] = obj.dataCodingInfo();
            %
            % Output Arguments:
            %
            % A --> lower limit for data {obj.Xlo}
            % B --> upper limit for data {obj.Xhi}
            % M --> median of data range
            %--------------------------------------------------------------
            A = obj.Xlo;
            B = obj.Xhi;
            M = mean( [ B; A ] );
        end % dataCodingInfo
    end % protected method signatures

    methods ( Access = private )
        function obj = setXlo( obj, A )
            %--------------------------------------------------------------
            % Set the low limit for the data coding
            %
            % obj = obj.setXlo( A );
            %
            % Input Arguments:
            %
            % A     --> (double) Lower limit for data coding: A --> [-1]
            %--------------------------------------------------------------
            if ( nargin < 2 ) || isempty( A ) || ( ( numel( A ) ~= obj.N )...
                                              && ~isempty( obj.X ) )
                %----------------------------------------------------------
                % Apply default
                %----------------------------------------------------------
                A = min( obj.X );
            end
            obj.Xlo = A;
        end % setXlo

        function obj = setXhi( obj, B )
            %--------------------------------------------------------------
            % Set the high limit for the data coding
            %
            % obj = obj.setXlo( B );
            %
            % Input Arguments:
            %
            % B     --> (double) Lower limit for data coding: B --> [+1]
            %--------------------------------------------------------------
            if ( nargin < 2 ) || isempty( B ) || ( ( numel( B ) ~= obj.N )...
                                              && ~isempty( obj.X ) )
                %----------------------------------------------------------
                % Apply default
                %----------------------------------------------------------
                B = max( obj.X );
            end
            obj.Xhi = B;
        end % setXhi
    end % private method signatures

    methods
        function F = get.Fmax( obj )
            % Return maximum function query value to date
            F = max( obj.Y );
        end

        function X = get.Xmax( obj )
            % Return location of best known function query to date
            [ ~, Idx ] = max( obj.Y );
            X = obj.X( Idx, : );
        end
        
        function F = get.Fmin( obj )
            % Return minimum function query value to date
            F = min( obj.Y );
        end

        function X = get.Xmin( obj )
            % Return location of best known function query to date
            [ ~, Idx ] = min( obj.Y );
            X = obj.X( Idx, : );
        end
        
        function N = get.N( obj )
            % Return number of predictor variables
            N = size( obj.X, 2 );
        end % get.N

        function N = get.NumPoints( obj )
            % Return number of data points in the training set
            N = size( obj.X, 1 );
        end % get.NumPoints

        function Ok = get.DataOk( obj )
            % Check consistency of data
            Ok = size( obj.X, 1 ) & numel( obj.Y );
            Ok = Ok & ~isempty( obj.X ) & ~isempty( obj.Y );
        end % get.DataOk

        function Xc = get.Xc( obj )
            % Return coded input data [a,b] --> [-1,1]
            Xc = obj.code( obj.X );
        end % get.Xc
    end

    methods ( Access = protected, Static = true )
    end % private and static method signatures
end % classdef
