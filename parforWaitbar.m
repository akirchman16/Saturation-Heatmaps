function parforWaitbar(waitbarHandle,iterations)
    persistent count h N
    if nargin == 2
        %Initializes the waitbar
        count = 0;
        h = waitbarHandle;
        N = iterations;
    else
        %updates waitbar
        if isvalid(h)   %checks that handle is valid or if it refers to a deleted object
            count = count+1;
            waitbar(count/N,h);
        end
    end
end