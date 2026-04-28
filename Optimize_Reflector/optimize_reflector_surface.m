function [best_surf, best_fitness, fitness_history] = ...
    optimize_reflector_surface(init_surf, params)
%==========================================================================
% optimize_reflector_surface  —— Reflector surface evolutionary optimization function
%
% Function: Optimize the reflector surface perturbation values using a
%           Differential Evolution (DE) algorithm to maximize the
%           correlation coefficient between the H3-plane magnetic field
%           distribution and the target field pattern.
%
% Algorithm flow (per generation):
%   1. Call calc_H_field for each individual to compute the H3 field
%   2. Call eval_field_correlation to compute fitness (Re(corr_coeff) x 100)
%   3. Sort by fitness; elite individual is carried over directly
%   4. Roulette-wheel selection of parents
%   5. DE crossover: child = best + Fe * (r1 - r2)
%   6. Random mutation: randomly select mut_num nodes and add random perturbation
%
% Input:
%   init_surf  — 1x35 real vector, initial reflector control node
%                perturbation values (unit: m)
%   params     — struct, optimization hyperparameters (see field descriptions below)
%
% params fields (all have default values):
%   params.n_generations  — number of generations (default 50)
%   params.pop_size       — population size (default 30)
%   params.surf_depth     — perturbation depth limit (default 1/3.8e-3 m)
%   params.Fe             — DE crossover weight (default 0.1)
%   params.Fm             — mutation strength (default 0.1)
%   params.mut_num        — number of nodes mutated per step (default 2)
%   params.fig_handle     — figure handle for convergence plot (passed from main, optional)
%
% Output:
%   best_surf       — 1x35, best reflector control nodes found
%   best_fitness    — scalar, final best fitness value
%   fitness_history — n_generations x 1 vector, best fitness per generation
%==========================================================================

%% ========================================================================
%  Part 1: Read parameters and set defaults
%% ========================================================================

n_generations = get_field(params, 'n_generations', 50);
pop_size      = get_field(params, 'pop_size',      30);
surf_depth    = get_field(params, 'surf_depth',    1 / 3.8 * 1e-3);
Fe            = get_field(params, 'Fe',            0.1);
Fm            = get_field(params, 'Fm',            0.1);
mut_num       = get_field(params, 'mut_num',       2);
fig_handle    = get_field(params, 'fig_handle',    []);

n_ctrl = numel(init_surf);   % number of control nodes = 35

%% ========================================================================
%  Part 2: Initialize population (all individuals start from the same initial surface)
%% ========================================================================

population = cell(pop_size, 1);
for i = 1 : pop_size
    population{i} = init_surf;
end

fitness_history = zeros(n_generations, 1);

%% ========================================================================
%  Part 3: Main evolutionary loop
%% ========================================================================

for gen = 1 : n_generations

    %----------------------------------------------------------------------
    % 3.1 Evaluate fitness of each individual
    %----------------------------------------------------------------------
    fitness_vector = zeros(pop_size, 1);
    for l = 1 : pop_size
        H3 = calc_H_field(population{l});
        fitness_vector(l) = eval_field_correlation(H3);
    end

    %----------------------------------------------------------------------
    % 3.2 Record the best fitness of this generation and update the convergence plot
    %----------------------------------------------------------------------
    fitness_history(gen) = max(fitness_vector);

    if ~isempty(fig_handle) && ishandle(fig_handle)
        figure(fig_handle);
        plot(1:gen, 100 - fitness_history(1:gen), 'b-o', ...
             'LineWidth', 1.5, 'MarkerSize', 4);
        xlabel('Generation');
        ylabel('Error (100 - Fitness)');
        title('Reflector Surface Optimization Convergence');
        grid on;
        drawnow;
    end

    fprintf('Gen %3d | Best fitness: %.4f\n', gen, fitness_history(gen));

    %----------------------------------------------------------------------
    % 3.3 Sort individuals by fitness in ascending order
    %     (last index corresponds to the best individual)
    %----------------------------------------------------------------------
    [~, sort_idx] = sort(fitness_vector, 'ascend');

    %----------------------------------------------------------------------
    % 3.4 Roulette-wheel selection: select parents proportional to fitness
    %----------------------------------------------------------------------
    % Rearrange population from highest to lowest fitness
    pop_desc     = cell(pop_size, 1);
    fitness_desc = zeros(pop_size, 1);
    for i = 1 : pop_size
        pop_desc{i}     = population{sort_idx(pop_size + 1 - i)};
        fitness_desc(i) = fitness_vector(sort_idx(pop_size + 1 - i));
    end

    % Normalize to selection probabilities and compute cumulative distribution
    select_prob = fitness_desc / sum(fitness_desc);
    cumul_prob  = cumsum(select_prob);

    % Elite individual is carried over directly (highest fitness)
    new_pop    = cell(pop_size, 1);
    new_pop{1} = pop_desc{1};

    % Roulette-wheel selection for remaining parents
    for k = 2 : pop_size
        r = rand();
        w = 1;
        while w <= pop_size
            if r <= cumul_prob(w)
                new_pop{k} = pop_desc{w};
                break;
            end
            w = w + 1;
        end
    end

    %----------------------------------------------------------------------
    % 3.5 DE crossover (DE/best/1 strategy)
    %  child = best + Fe * (r1 - r2)
    %  best is the current best individual; r1 and r2 are two randomly
    %  selected distinct non-elite individuals
    %----------------------------------------------------------------------
    children    = cell(pop_size, 1);
    children{1} = new_pop{1};   % elite individual passes directly to next generation

    best_individual = new_pop{1};

    for mm = 2 : pop_size
        i_r1 = round(rand() * (pop_size - 2)) + 2;
        i_r2 = round(rand() * (pop_size - 2)) + 2;

        child = best_individual + Fe * (new_pop{i_r1} - new_pop{i_r2});

        % Boundary constraints
        child = max(child, -surf_depth / 2);
        child = min(child,  surf_depth / 2);

        children{mm} = child;
    end

    %----------------------------------------------------------------------
    % 3.6 Random mutation
    %  For each non-elite individual, randomly select mut_num nodes and
    %  add a random perturbation scaled by Fm
    %----------------------------------------------------------------------
    for k = 2 : pop_size
        mut_indices = unidrnd(n_ctrl, [mut_num, 1]);
        mut_delta   = (rand(mut_num, 1) - 0.5) * surf_depth;

        child_m = children{k};
        for xy = 1 : mut_num
            new_val = child_m(mut_indices(xy)) + mut_delta(xy) * Fm;
            new_val = max(new_val, -surf_depth / 2);
            new_val = min(new_val,  surf_depth / 2);
            child_m(mut_indices(xy)) = new_val;
        end
        population{k} = child_m;
    end
    population{1} = children{1};   % preserve elite individual

end

%% ========================================================================
%  Part 4: Return the best result
%% ========================================================================

% Evaluate fitness one final time on the last generation to identify the best individual
fitness_final = zeros(pop_size, 1);
for l = 1 : pop_size
    H3 = calc_H_field(population{l});
    fitness_final(l) = eval_field_correlation(H3);
end

[best_fitness, best_idx] = max(fitness_final);
best_surf = population{best_idx};

end

%% ========================================================================
%  Helper function: safely read a struct field, returning a default value
%  if the field does not exist
%% ========================================================================

function val = get_field(s, field_name, default_val)
    if isfield(s, field_name)
        val = s.(field_name);
    else
        val = default_val;
    end
end
