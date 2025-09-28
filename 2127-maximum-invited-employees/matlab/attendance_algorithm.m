%% Notes

% TODO
% - edge case handling: 
%   - n = 1
%   - oneself as favorite

% optimization (optional)
% - previous valid loop detection, i.e. terminate iteration on loop detect
% - specify maximum attending case table(s) (not required)

%% Init

clear; clc;

%% Input

% data specification
method = "manual";  % "manual" or "generated"

% debugging
dbg = 0;

%% Test Cases

N = 18;

example_index = 1;

test_data = {
    
    % pronlem examples
    [2, 2, 1, 2];
    [1, 2, 0];
    [3, 0, 1, 4, 1];

    % custom examples
    [4, 4, 5, 2, 2, 1];
    [4, 4, 5, 6, 3, 2, 4, 8, 7];
    [1, 0, 1, 1, 3, 3, 3];
    [2, 2, 1, 2, 1, 4, 0];
    [2, 2, 1, 2, 1, 4, 0, 8, 7, 8, 9, 12, 11];

    % further provided examples
    [1, 2, 3, 4, 5, 0];
    [6, 4, 4, 5, 0, 3, 3];
    [1, 2, 3, 4, 5, 6, 3, 8, 9, 10, 11, 8];
    [1, 0, 0, 2, 1, 4, 7, 8, 9, 6, 7, 10, 8];
    [1, 0, 3, 2, 5, 6, 7, 4, 9, 8, 11, 10, 11, 12, 10];
    [44, 47, 86, 29, 65, 58, 0, 86, 45, 96, 21, 41, 59, 95, 55, 25, 58, 6, 19, 56, 66, 31, 76, 50, 72, 14, 59, 81, 59, 83, 97, 19, 36, 41, 32, 6, 18, 30, 56, 59, 77, 24, 35, 60, 79, 29, 26, 26, 49, 15, 93, 53, 9, 30, 96, 20, 51, 72, 53, 87, 85, 82, 53, 8, 24, 56, 37, 33, 5, 96, 83, 61, 65, 67, 54, 88, 69, 81, 82, 55, 10, 19, 40, 56, 73, 6, 72, 28, 65, 4, 40, 94, 49, 57, 34, 71, 85, 35];
    [131, 316, 117, 308, 345, 374, 103, 135, 319, 422, 376, 58, 175, 269, 17, 265, 375, 374, 419, 77, 404, 190, 48, 2, 25, 91, 37, 150, 42, 119, 182, 62, 408, 24, 359, 358, 17, 243, 243, 53, 271, 48, 395, 45, 367, 38, 96, 44, 302, 143, 351, 360, 163, 70, 353, 391, 224, 201, 366, 381, 173, 414, 92, 57, 398, 36, 338, 219, 418, 221, 324, 314, 335, 353, 88, 120, 54, 117, 379, 200, 240, 92, 83, 306, 81, 236, 420, 9, 195, 70, 152, 156, 162, 410, 288, 344, 125, 9, 217, 359, 84, 161, 195, 377, 174, 136, 22, 351, 71, 411, 319, 77, 228, 368, 264, 80, 156, 47, 365, 34, 325, 133, 174, 169, 237, 140, 176, 421, 21, 72, 301, 210, 286, 212, 316, 247, 287, 93, 350, 117, 372, 374, 283, 192, 130, 0, 276, 265, 122, 266, 420, 416, 360, 83, 21, 126, 250, 97, 246, 382, 376, 388, 378, 286, 403, 256, 322, 411, 238, 211, 25, 345, 361, 322, 222, 37, 417, 363, 258, 286, 170, 118, 148, 272, 151, 83, 213, 127, 207, 338, 7, 277, 124, 238, 8, 147, 190, 117, 231, 60, 6, 254, 187, 83, 193, 239, 312, 325, 287, 404, 191, 389, 234, 152, 149, 45, 397, 136, 74, 77, 265, 216, 389, 330, 90, 380, 385, 358, 15, 379, 339, 162, 13, 285, 120, 30, 413, 49, 328, 181, 356, 259, 149, 316, 163, 61, 301, 85, 118, 152, 217, 21, 118, 244, 345, 363, 104, 240, 3, 46, 128, 219, 410, 185, 186, 96, 385, 245, 4, 79, 228, 362, 346, 309, 90, 228, 132, 94, 280, 35, 241, 207, 305, 229, 53, 154, 212, 336, 307, 397, 59, 99, 152, 115, 408, 333, 21, 176, 320, 24, 212, 328, 96, 237, 169, 368, 356, 141, 386, 219, 25, 189, 370, 232, 334, 181, 372, 139, 277, 5, 348, 251, 15, 28, 145, 187, 51, 83, 112, 300, 300, 215, 104, 217, 393, 215, 242, 16, 134, 352, 420, 45, 291, 132, 151, 422, 369, 335, 114, 318, 362, 287, 192, 76, 39, 52, 305, 284, 307, 223, 326, 74, 296, 212, 99, 101, 176, 411, 50, 16, 392, 421, 8, 322, 143, 330, 122, 389, 259, 59, 226, 156, 372, 393, 349, 223, 357, 270, 182, 114, 82, 318, 75, 182, 123, 38, 236, 296, 41, 14, 293, 100, 320, 87, 158, 276, 353, 77, 103, 50, 321, 422, 274, 299, 250, 183, 45, 374, 24, 124, 260, 371, 306, 217, 75];

    % special cases
    [1:N, 0];  % all attending 
    [1, 0];  % pair-only
    1;  % "n = 1"

};

%% Data

if method == "manual"  % manual

    favorites = test_data{example_index};
    
elseif method == "generated"  % generated

    s = rng;

    indices = (0:(N - 1));
    list = randi(N, 1, N) - 1;
    idx = list(list == indices);

    for i = 1:length(idx)

        sd = setdiff(0:(N - 1), idx(i));
        list(idx(i) + 1) = sd(randi((N - 1) - 1));

    end

    favorites = list;

end

%% Algorithm

% data
n = numel(favorites);
employees = (0:(n - 1));
mapping = [employees; favorites];

if dbg == 1

    mapping

end

% pre-allocation
seq = zeros(n);  % working sequence matrix

seq_lengths = zeros(n, 1);
pair_roots = zeros(n, 2);  % pair-loops (potential sequence roots)

seq_length = 0;
pair_only_length = 0;
max_loop_length = 0;

% independent indices
q = 1; p = 1;

for i = 1:n  % MATLAB indexing (problem: 0-indexed)

    if n < 2

        disp("Invalid Input Data")

        return  % terminate script exection

    end

    seq_case = zeros(1, n);  % table array

    j = i;  % indexing

    cl = 1;  % column index

    while true

        e = employees(j);
        
        if j == i

            seq_case(cl) = e;
            cl = cl + 1;

        end      

        f = favorites(e + 1);  % shifted index -> 1-indexed (MATLAB)

        if ismember(f, seq_case(1:(cl - 1)))

            break

        else

            seq_case(cl) = f;
            cl = cl + 1;

        end

        j = f + 1;

    end
    
    seq_length = cl - 1;
    
    %%% Validity Check:
    %
    % 1. "Valid Loop": valid table with closed loop attendance
    %   - last empl.'s favorite is first empl. on table (incl. singular pairs)
    %       -> logic: t(1) == favorites( t(end) + 1 )
    %
    % 2. Table larger than single pair ends with pair (Pair-stall)
    %   - last empl.'s favorite is preceding sequence entry
    %       -> logic: t(end - 1) == favorites( t(end) + 1 )
    %%%

    % "Valid Loop"
    if ( favorites( seq_case(seq_length) + 1 ) == seq_case(1) )  

        if seq_length > max_loop_length

            max_loop_length = seq_length;

        end

        if seq_length == 2

            pair_only_length = pair_only_length + 1;
            pair_roots(p, :) = seq_case(1:seq_length);

            p = p + 1;

        else

            pair_roots(end, :) = [];

        end

    % efficiency
    seq(end, :) = [];
    seq_lengths(end) = [];

    % Pair Stall
    elseif ( favorites( seq_case(seq_length) + 1 ) == seq_case(seq_length - 1) )

        seq(q, :) = seq_case(:);
        seq_lengths(q) = seq_length;

        q = q + 1;

        % efficiency
        pair_roots(end, :) = [];

    else

        % efficiency
        seq(end, :) = [];
        seq_lengths(end) = [];
        pair_roots(end, :) = [];

    end
    
end

max_pair_only_length = pair_only_length;

if dbg == 1

    seq
    seq_lengths
    pair_roots
    max_loop_length
    max_pair_only_length

end

for ii = 1:size(pair_roots, 1)

    qq = 1;

    tmp = [];

    for jj = 1:size(seq, 1)

        sl = seq_lengths(qq);

        % check root pair
        if isequal( seq(qq, (sl - 1):sl), pair_roots(ii, :) )

            if (length(seq(qq, 1:sl)) >= length(tmp)) && (~isempty(tmp))

                    % omit processed data (efficiency)
                    seq(pp, :) = [];
                    seq_lengths(pp) = [];

                    continue

            end

            pp = qq;
            tmp = seq(pp, 1:sl);

        end

        qq = qq + 1;

    end

end

if dbg == 1

    seq

end

max_pair_root_chain_length = 0;

if ~isempty(seq_lengths)

    max_pair_root_chain_length = max_pair_only_length + sum(seq_lengths - 2);

    if dbg == 1

        max_pair_root_chain_length

    end

end

%% Output

table_length_cases = [

        max_loop_length;
        max_pair_only_length;
        max_pair_root_chain_length;
        
];

if dbg == 1

    table_length_cases

end

max_attendance = max(table_length_cases);

disp("Max. Attending Employees: " + num2str(max_attendance))