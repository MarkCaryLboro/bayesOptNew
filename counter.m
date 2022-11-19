classdef counter < handle
    %----------------------------------------------------------------------
    % A class to implement a mod-based counter
    %----------------------------------------------------------------------

    properties ( SetAccess = protected )
        Counter     (1,1)   int64   { mustBeGreaterThanOrEqual( Counter, 0 ) } = 0
        Base        (1,1)   int64   { mustBeGreaterThan( Base, 1 ) } = 2
        Ratchet     (1,1)   int64   { mustBeGreaterThanOrEqual( Ratchet, 0 ) } = 0
    end % protected properties

    properties ( Access = private, Dependent = true )
        Counter_ 
        Base_
        Ratchet_
    end % private properties

    methods
        function obj = counter( Initial, Base, Ratchet)
        end % class constructor
    end % ordinary methods signatures

    methods
        function C = get.Counter_( obj )
            C = double( obj.Counter );
        end

        function B = get.Base_( obj )
            B = double( obj.Base );
        end

        function R = get.Ratchet_( obj )
            R = double( obj.Ratchet );
        end
    end % Set/Get methods
end % classdef