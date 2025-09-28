%% Notes

% optimize?
% edge cases: n = 1, oneself as favorite
% - OPTIONAL: specify table(s) of maximum attending case (not required)

%% Init

clear; clc;

s = rng;

%% Input

% Test Cases

% 1. All attending
% N = 18;
% test_data = (0:N) + 1;
% test_data(end) = 0;

% 2. Problem examples
% test_data = [2, 2, 1, 2];  % Example 1
% test_data = [1, 2, 0];  % Example 2
% test_data = [3, 0, 1, 4, 1];  % Example 3

% 3. Custom examples
% test_data = [4, 4, 5, 2, 2, 1];
 test_data = [4, 4, 5, 6, 3, 2, 4, 8, 7];
% test_data = [1, 0, 1, 1, 3, 3, 3];
% test_data = [2, 2, 1, 2, 1, 4, 0];
% test_data = [2, 2, 1, 2, 1, 4, 0, 8, 7, 8, 9, 12, 11]; 

% Data specification
method = "man";  % 'man' (manual) or 'gen' (generated)

%% Data

if method == "man"  % manual

    favorites = test_data;  % note: 0-indexed notation specified
    n = numel(favorites);
    
elseif  method == "gen"  % generated

    n = 19;

    indices = (0:(n - 1));
    list = randi(n, 1, n) - 1;
    idx = list(list == indices);

    for i = 1:length(idx)

        sd = setdiff(0:(n - 1), idx(i));
        list(idx(i) + 1) = sd(randi((n - 1) - 1));

    end

    favorites = list;

end

%% Algorithm

employees = (0:(n - 1));
mapping = [employees; favorites];
mapping

seq = zeros(n);  % sequence matrix

seq_length = 0;
seq_lengths = zeros(n, 1);
pair_roots = zeros(n, 2);

pair_only_length = 0;
max_loop_length = 0;

q = 1;
p = 1;

for i = 1:n  % MATLAB notation; default: 0-indexed problem: 0:(n-1)

    if n < 2

        disp("Invalid Input Data.")
        break

    end

    seq_case = zeros(1, n);  % circular table array

    j = i;  % indexing

    cl = 1;  % column iterator

    while true

        e = employees(j);
        
        if j == i  % or cl == 1

            seq_case(cl) = e;
            cl = cl + 1;

        end      

        f = favorites(e + 1);  % shift -> 1-indexed (MATLAB)

        if ismember(f, seq_case(1:(cl - 1)))

            break

        else

            seq_case(cl) = f;
            cl = cl + 1;

        end

        j = f + 1;

    end

    % disp("--------------------------------------------")
    % seq_case
    
    seq_length = cl - 1;
    
    % Check row validity:
    %
    % 1. Valid table loop (Loop) 
    % Last empl.'s favorite is first empl. on table (includes single pairs)
    %       -> logic: t(end) == favorites(t(1) + 1) )
    %
    % 2. Table larger than single pair ends with pair (Pair-stall)
    % Last empl.'s favorite is ONLY preceding entry
    %       -> logic: t(end) == favorites(t(end - 1) + 1))
    %
    % NOTE: `+ 1` due to MATLABs 1-indexed arrays (0th employee present)

    % WRONG if ( seq_case(seq_length) == favorites(seq_case(1) + 1) )  % Loop (valid)
    if ( favorites(seq_case(seq_length) + 1) == seq_case(1) )  % Loop (valid)

        % disp("Valid Loop")
        % continue

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

    seq(end, :) = [];
    seq_lengths(end) = [];

    % WRONG elseif ( seq_case(seq_length) == favorites(seq_case(seq_length - 1) + 1) )  % Pair-stall
    elseif ( favorites(seq_case(seq_length) + 1) == seq_case(seq_length - 1) )  % Pair-stall

        % disp("Pair Stall")
        % continue

        seq(q, :) = seq_case(:);
        seq_lengths(q) = seq_length;

        q = q + 1;

        pair_roots(end, :) = [];

    else

        seq(end, :) = [];
        seq_lengths(end) = [];
        pair_roots(end, :) = [];

    end
    
end


% now only contain non-loop non-pair-only valid cases
seq
seq_lengths  
pair_roots

max_loop_length
max_pair_only_length = pair_only_length

for ii = 1:size(pair_roots, 1)

    qq = 1;

    tmp = [];

    for jj = 1:size(seq, 1)

        sl = seq_lengths(qq);

        % check root pair
        if isequal( seq(qq, (sl - 1):sl), pair_roots(ii, :) )

            if (length(seq(qq, 1:sl)) >= length(tmp)) && (~isempty(tmp))

                    % eliminate processed data for efficiency
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

seq

max_pair_root_chain_length = 0;
if ~isempty(seq_lengths)
    max_pair_root_chain_length = max_pair_only_length + sum(seq_lengths - 2);
    max_pair_root_chain_length
end


%% Output

table_length_cases = [

        max_loop_length,
        max_pair_only_length,
        max_pair_root_chain_length

    ]

max_attendance = max(table_length_cases);

disp("Maximum number of attending employees: " + num2str(max_attendance))