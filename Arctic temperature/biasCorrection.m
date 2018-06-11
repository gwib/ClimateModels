function [v1_corrected,std] = biasCorrection(v1)
    std = NaN(320,160);
    for i = 1:320
        for j=1:160
            y_v1 = v1(i,j,1:10);
            %disp(y_v1)
            %disp(max(max(y_v1)))
            delta_T = max(max(y_v1)) - min(min(y_v1));
            std(i,j) = (nanmean(y_v1) - min(min(y_v1)))/delta_T;
        end
    end
    v1_corrected = v1 - std;
end