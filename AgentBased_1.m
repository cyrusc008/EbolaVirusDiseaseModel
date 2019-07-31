%% Agent-based Model V.1 (CYRUS)

clear all, clc

%initialize the lattice
lattice = zeros(10,10);

%seed agents (cells)
n_init = 10;
order = randperm(100);
sites = order(1:20);
lattices(sites) = 1;

%thresholds
p_division = 0.2;
p_movement = 0.1;

%begin simulation and initialize storage
t = 1;
t_limit = 10;
lattice_storage = NaN(10,10,t_limit);
lattice_storage(:,:,1) = lattice;

while t+1 <= t_limit
    for i = 1:100
        
        %proceed with rules if agent is located in site i
        if lattice(i) == 1
            %get neighborhood around current agent 
            current_agent_position = zeros(10,10);
            current_agent_position(i) = 1;
            distance_matrix = bwdist(current_agent_position,'chessboard');
            
            %rule 1: if there is at least one empty neighboring spot, the
            %cell may divide with probability p_division 
            f = find(distance_matrix == 1);
            empty_spots = f(lattice(f) == 0);
            
            if numel(empty_spots) >= 1
                idx = randi(numel(empty_spots));
                q = rand(1);
                if q <= p_division
                    lattice (empty_spots(idx)) = 1;
                end
            end
            
            %rule 2: if there is at least one empty neighboring spot, the
            %cell may move with probability p_movement
            empty_spots = f(lattice(f) == 0);
            if numel(empty_spots) >= 1
                idx = randi(numel(empty_spots));
                q = rand(1);
                if q <= p_movement
                    lattice(empty_spots(idx)) = 1;
                    lattice(i) = 0;
                end
            end
        end
    end
    
    %store current lattice and update count
    lattice_storage(:,:,t+1) = lattice;
    t = t+1;
end