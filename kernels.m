classdef kernels < int8
    % define permitted kernel type for GPR

    enumeration
        ARDsquaredExponential    (0)
        ARDexponential           (1)
        ARDmatern32              (2)
        ARDmatern52              (3)
    end
end % classdef