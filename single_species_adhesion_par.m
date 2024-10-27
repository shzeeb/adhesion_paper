%adhesion_v1_par.m: single species model with cell-cell adhesion on a two
%dimensional lattice

%Coded by: Shahzeb Raja Noureen

clear;

%number of lattice points in the x (ncols) and y (nrows) direction
nrows=20;
ncols=400;

%movement rate
r_m=1;

%adhesion strength
q=0.25;

%simulation run time
T_final=1000;

%Initial condition
IC=2;%set IC=2. See initial_cond() function for more details.

%number of repetitions
reps=50;

%recording step
rec_step=100;

%number of recording steps
n_rec_steps = T_final/rec_step;

%recording matrix for averaging cell density
rec_mat_full=zeros(nrows,ncols,n_rec_steps+1,reps);

%recording times for video
rec_times = 1:rec_step:T_final;

len_rt = length(rec_times);


%replace for with parfor for parallel computing
parfor i=1:reps

%initialisng the domain
A=initial_cond(nrows,ncols,IC);

%recording matrix
rec_mat = zeros(nrows,ncols,n_rec_steps+1);

rec_mat(:,:,1)=A;
%recoding matrix A for movie
% A_rec_mat=zeros(nrows,ncols,len_rt);

%initial domain matrix
A_init = A;

%recording times for video
rec_times = 1:1:T_final;

len_rt = length(rec_times);


%initial state of the recording matrix
% rec_mat(:,:,1) = A_init;

t_after = 0;

%find the initial agent positions for both species
[x,y]=find(A==1);


%number of M and X species
M=length(x);


%propensities
a=r_m*M;

%sum of propensities
a0 = a;

%initialise time
t = 0;

%recording index (for movie)
j = 2;

k=0;

% the main loop
while t < T_final
    k=k+1;
    
    r=rand;
    
    %generate time to next event
    tau = log(1/r)/a0;
    
    %update current time
    t = t + tau;   
    
    %generate uniformly disttributed random number
    r = rand;
    
    ra0=r*a0;

    z=[x,y];
    
    %cell type for movement
    ind_agent=randi(size(z,1));
    agent=z(ind_agent,:);

    %%need to reimplement this
    %neighbouring sites
    neighbours = [agent(1)-1 agent(2); agent(1)+1 agent(2); agent(1) agent(2)-1; agent(1) agent(2)+1];
    
    %draw a target size at random with probability 1/4
    target_site = datasample(neighbours,1,1);
    
    %check for boundaries
    if isempty(target_site) || target_site(1) == 0 || target_site(2) == 0 || target_site(1) == nrows+1 || target_site(2) == ncols+1
%         disp('no move');
    else

        %check movement if target_site empty and cell_type==1
        if A(target_site(1),target_site(2)) == 0
            
            %find sum of neighbour
            sum_neigh=0;
            for kk=1:4
                if neighbours(kk,1) ~= 0 && neighbours(kk,2) ~= 0 && neighbours(kk,1) ~= nrows+1 && neighbours(kk,2) ~= ncols+1
                    sum_neigh=sum_neigh+A(neighbours(kk,1),neighbours(kk,2));
                end
            end
            
            if sum_neigh==0
                A(target_site(1),target_site(2)) = A(agent(1),agent(2));
                A(agent(1),agent(2)) = 0;

                x(ind_agent)=target_site(1);
                y(ind_agent)=target_site(2);
                
            else
            
                %movement probabilitty due to adhesion
                move_prob=(1-q)^sum_neigh;

                r=rand;

                %choose to move/not move based on movement probability
                if r<=move_prob

                    A(target_site(1),target_site(2)) = A(agent(1),agent(2));
                    A(agent(1),agent(2)) = 0;

                    x(ind_agent)=target_site(1);
                    y(ind_agent)=target_site(2);
                end
            end
        end
    end
    
    
%calculate the times for recording
t_before = t_after;
t_after = t;

ind_before = ceil((t_before+eps)/rec_step);
ind_after = min(floor(t_after/rec_step),n_rec_steps);

%steps to write to:
steps_to_write = ind_after-ind_before+1;

if steps_to_write > 0 && steps_to_write ~= inf

    if steps_to_write > 1

        rec_mat(:,:,ind_before+1:ind_after+1)=rec_mat(:,:,ind_before+1:ind_after+1)+repmat(A,1,1,size(rec_mat(:,:,ind_before+1:ind_after+1),3));

    else
        rec_mat(:,:,ind_before+1:ind_after+1)=rec_mat(:,:,ind_before+1:ind_after+1)+A;
    end
end

end

rec_mat_full(:,:,:,i)=rec_mat;
end


%save workspace
save("adhesion_IC="+num2str(IC)+"_q="+num2str(q)+"_T="+num2str(T_final)+".mat",'-v7.3')