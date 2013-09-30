function [ lMin, lMax ] = gershgorin (A)
    lMin = A(1, 1);
    lMax = A(1, 1);
    for i = 1:size(A, 1)
        R = 0;
        for j = 1:size(A, 2)
            if(i ~= j)
                R = R+abs(A(i, j));
            end
        end
        if(A(i, i)-R < lMin)
            lMin = A(i, i)-R;
        end
        if(A(i, i)+R > lMax)
            lMax = A(i, i)+R;
        end
    end
