function drawPlot( B, lims, spec )
%drawPlot draw line y=b(1)+x*B(2) in axis with limits 'lims' by line
%specification 'spec'.
    
    % Calculate value of Y i n the two borders of x interval
    x1 = lims(1);
    x2 = lims(2);
    y1 = B(1) + x1 * B(2);
    y2 = B(1) + x2 * B(2);
    if y1 > y2
        % Decreasing line
        if y1 > lims(4)
            y1 = lims(4);
            x1 = (y1 - B(1)) / B(2);
        end
        if y2 < lims(3)
            y2 = lims(3);
            x2 = (y2 - B(1)) / B(2);
        end
    else
        % Increasing line
        if y2 > lims(4)
            y2 = lims(4);
            x2 = (y2 - B(1)) / B(2);
        end
        if y1 < lims(3)
            y1 = lims(3);
            x1 = (y1 - B(1)) / B(2);
        end
    end
    % Draw line
    plot([x1, x2], [y1, y2], spec);
end

