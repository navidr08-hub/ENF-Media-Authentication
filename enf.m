function [y1, y2, y3] = enf(x, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency)
    % Validate inputs
    expectedSize = [NaN, 1];
    validateattributes(x, {'double'}, {'vector', 'nonempty', 'size', expectedSize});
    validateattributes(Fs, {'numeric'}, {'scalar', 'integer'});
    validateattributes(BlockSize, {'numeric'}, {'scalar', 'integer'});
    validateattributes(ZeroPad, {'numeric'}, {'scalar', 'integer'});
    validateattributes(Window, {'double'}, {'vector', 'nonempty', 'size', expectedSize});

    if ~isscalar(Frequency) && ~isvector(Frequency)
        error('Input must be a scalar or a vector.');
    else
        validateattributes(Frequency, {'numeric'}, {'integer'});
    end

    % Segmentation
    increment =(BlockSize - (BlockSize*Overlap));
    data_size = (length(x) - BlockSize);
    incr_indeces = (0:BlockSize-1);
    cols = (1:increment:data_size+1);
    IndexMatrix = repmat(cols, [BlockSize, 1]) + repmat(incr_indeces', [1, length(cols)]);
    Segmented = x(IndexMatrix);
    
    % Windowing
    Hann_Window = repmat(Window, [1, length(cols)]);
    Windowed = Segmented.*Hann_Window;
    
    % Zero-Padding
    Padded = padarray(Windowed, ZeroPad, 'post');
    
    % Fast-Fourier Transform
    N = BlockSize + ZeroPad;
    FastFT = abs(fft(Padded,N,1));

    % Frequency Extraction
    function [range] = desired_ranges(f, Frequency)
        [~, c] = size(Frequency);
        if (c > 1)
            range = zeros(1, length(f));
            for i = Frequency
                range = range + (f>=(i-1) & f<=(i+1));
            end
        else
            range = f>=(Frequency-1) & f<=(Frequency+1);
        end
    end

    f = Fs*[0:N-1]/N;
    k = find(desired_ranges(f, Frequency));
    Fk = f(k);
    
    % Magnitude Plot
    y1 = abs(FastFT(k, :));                         % need to transpose
    
    % Maximum Magnitude Calculations
    [~,max_mag] = max(y1,[],1);
    y2 = Fk(max_mag);

    % Weighted Magnitude Calculations
    desired_freq = repmat(Fk', 1, length(cols));
    numerator = sum(desired_freq.*y1, 1);
    denominator = sum(y1, 1);
    y3 = numerator ./ denominator;                  % need to transpose
end