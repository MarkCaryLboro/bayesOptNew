classdef counter < handle
    %----------------------------------------------------------------------
    % A class to implement a mod-based counter
    %----------------------------------------------------------------------

    events
    end

    properties ( SetAccess = protected )
        Counter     (1,1)   int64   { mustBeGreaterThanOrEqual( Counter, 0 ) } = 0
        Base        (1,1)   int64   { mustBeGreaterThan( Base, 1 ) } = 2
        Ratchet     (1,1)   int64   { mustBeGreaterThanOrEqual( Ratchet, 0 ) } = 1
    end % protected properties

    properties ( Access = private, Dependent = true )
        Counter_ 
        Base_
        Ratchet_
    end % private properties

    methods
        function obj = counter( Initial, Base, Ratchet)
            %--------------------------------------------------------------
            % Class constructor. Define a counter object. The counter is
            % updated when rem( Counter, Base ) == Ratchet.
            %
            % obj = counter( Initial, Base, Ratchet);
            %
            % Input Arguments:
            %
            % Initial   --> Initial counter value {0}
            % Base      --> Number divisor
            % Ratchet   --> Increment the counter when the remainder is
            %               equal to this value.
            %--------------------------------------------------------------
            arguments
                Initial     (1,1)   int64   { mustBeGreaterThanOrEqual( Initial, 0 ) } = 0
                Base        (1,1)   int64   { mustBeGreaterThan( Base, 1 ) } = 2
                Ratchet     (1,1)   int64   { mustBeGreaterThanOrEqual( Ratchet, 0 ) } = 1
            end
            obj.Counter = Initial;
            obj.Base = Base;
            obj.Ratchet = Ratchet;
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