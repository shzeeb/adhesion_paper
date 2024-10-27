%two species differential adhesion in two dimensions. ABM stochastic
%simulation based on the gillespie algorithm.

%Coded by: Shahzeb Raja Noureen


clear;

%number of lattice points in the x (ncols) and y (nrows) direction
nrows=20;
ncols=200;

%movement rate type M
r_a=1;

%movement rate type X
r_b=1;

%swapping probability
rho=1;

%video indicator: set to 1 --> record video; 0 --> do not record video
video=0;

%save data indicator: set to 1 --> save data; 0 --> do not save data
Save=0;

%adhesion strength type A
p=0;

%adhesion strength type B
q=0;

%cross-adhesion strength
r=0;

T_final=1000;

%type of initial condition. See function called initial_cond() to choose
%the appropriate IC
IC=1;

%number of repetitions
nreps=10;

%recording step
rec_step=100;

%recording times
rec_times = 0:rec_step:T_final;

%number of recording steps
len_rt = length(rec_times);

%for video
if video==1
    v = VideoWriter('adhesion_2D_video.avi');
    v.FrameRate = 20;
    v.Quality = 100;
    open(v);
    
    set(gcf,'position',[0,0,1080,1920]);
    
    c=zeros(3,3);
    
    c(1,:)=[1 1 1];
    c(2,:)=[1 0 0];
    c(3,:)=[0 166/255 81/255];
    
end

%recording matrix for averaging cell density
rec_mat_full=zeros(nrows,ncols,len_rt,nreps);

parfor i=1:nreps
    
    rng('shuffle');
    
    
    %initialisng the domain
    A=initial_cond(nrows,ncols,IC);
    
    
    %recording times for video
    rec_times=0:rec_step:T_final;
    
    %recording matrix
    rec_mat=zeros(nrows,ncols,len_rt);
    
    %initial domain matrix
    A_init = A;
    
    rec_mat(:,:,1)=A_init;
    
    %initialise time
    t = 0;
    
    %recording index (for movie)
    j = 2;
    
    %number of M and X species
    N_a=sum(sum(A==1));
    N_b=sum(sum(A==2));
    
    % the main loop
    while t < T_final
        
        %find agent indices
        [x1,y1]=find(A==1);
        [x2,y2]=find(A==2);
        
        z1=[x1,y1];
        z2=[x2,y2];
        
        %propensities
        a=zeros(1,2);
        
        a(1)=r_a*N_a; %movement species A
        a(2)=r_b*N_b; %movement species B
        
        cumsuma=cumsum(a);
        
        a0=cumsuma(end);
        
        rnd=rand;
        
        %generate time to next event
        tau = -log(rnd)/a0;
        
        %update current time
        t = t + tau;
        
        rnd=rand;
        
        rnda0=rnd*a0;
        
        %sample agent and decide the type of event of occur
        if rnda0<=cumsuma(1)
            ind_agent=randi(size(z1,1));
            agent=z1(ind_agent,:);
            
        else
            ind_agent=randi(size(z2,1));
            agent=z2(ind_agent,:);
            
        end
        
        %neighbouring sites
        neighbours_agent = [agent(1)-1 agent(2); agent(1)+1 agent(2); agent(1) agent(2)-1; agent(1) agent(2)+1];
        
        %draw a target size at random with probability 1/4
        target_site = datasample(neighbours_agent,1,1);
        
        %check for boundaries
        if isempty(target_site) || target_site(1) == 0 || target_site(2) == 0 || target_site(1) == nrows+1 || target_site(2) == ncols+1
            
        else
                
            %check movement if target_site empty and cell_type==1
            if A(target_site(1),target_site(2))==0 && A(agent(1),agent(2))==1
                
                %find sum of neighbour
                sum_neigh_agent_A=0;
                sum_neigh_agent_B=0;
                
                for kk=1:4
                    if neighbours_agent(kk,1) ~= 0 && neighbours_agent(kk,2) ~= 0 && neighbours_agent(kk,1) ~= nrows+1 && neighbours_agent(kk,2) ~= ncols+1
                        if A(neighbours_agent(kk,1),neighbours_agent(kk,2))==1
                            sum_neigh_agent_A=sum_neigh_agent_A+1;
                            
                        elseif A(neighbours_agent(kk,1),neighbours_agent(kk,2))==2
                            sum_neigh_agent_B=sum_neigh_agent_B+1;
                            
                        end
                    end
                end
                
                move_prob_A=(1-p)^sum_neigh_agent_A;
                move_prob_B=(1-r)^sum_neigh_agent_B;
                
                move_prob=move_prob_A*move_prob_B;
                
                r1=rand;
                
                if r1<=move_prob
                    A(target_site(1),target_site(2))=1;
                    A(agent(1),agent(2))=0;
                    
                end
                
            elseif A(target_site(1),target_site(2)) == 0 && A(agent(1),agent(2))==2
                %if target site empty and the focal agent is a type B
                %agent
                
                %find sum of neighbour
                sum_neigh_agent_A=0;
                sum_neigh_agent_B=0;
                for kk=1:4
                    if neighbours_agent(kk,1) ~= 0 && neighbours_agent(kk,2) ~= 0 && neighbours_agent(kk,1) ~= nrows+1 && neighbours_agent(kk,2) ~= ncols+1
                        if A(neighbours_agent(kk,1),neighbours_agent(kk,2))==1
                            sum_neigh_agent_A=sum_neigh_agent_A+1;
                            
                        elseif A(neighbours_agent(kk,1),neighbours_agent(kk,2))==2
                            sum_neigh_agent_B=sum_neigh_agent_B+1;
                        end
                    end
                end
                
                %movement probabilities of type A and type B agent
                move_prob_B=(1-q)^sum_neigh_agent_B;
                move_prob_A=(1-r)^sum_neigh_agent_A;
                
                %overall movement probability
                move_prob=move_prob_A*move_prob_B;
                
                
                r1=rand;
                
                if r1<=move_prob
                    A(target_site(1),target_site(2))=2;
                    A(agent(1),agent(2))=0;
                    
                end
                
                %check for swapping if target_site already occupied (only
                %check for swap with agent at target site different from
                %agent at the focal site)
            elseif A(target_site(1),target_site(2))~=0 && rho>0 && A(agent(1),agent(2))~=A(target_site(1),target_site(2))
                
                %generate random number
                rnd=rand;
                
                %check for a swap for type A agent
                if rnd<=rho && A(agent(1),agent(2))==1
                    
                    %find sum of neighbour of the focal agent
                    sum_neigh_agent_A=0;
                    sum_neigh_agent_B=0;
                    
                    %find sum of neighbour of the target agent
                    sum_neigh_targ_A=0;
                    sum_neigh_targ_B=0;
                    
                    %target site neighbours
                    targ_site_neighs=[target_site(1)-1 target_site(2);target_site(1)+1 target_site(2);...
                        target_site(1) target_site(2)-1;target_site(1) target_site(2)+1];
                    
                    %loop through each target site neighbour and
                    %calculate the number of type A and type B
                    %neighbours
                    for kk=1:4
                        if targ_site_neighs(kk,1) ~= 0 && targ_site_neighs(kk,2) ~= 0 && targ_site_neighs(kk,1) ~= nrows+1 && targ_site_neighs(kk,2) ~= ncols+1
                            if A(targ_site_neighs(kk,1),targ_site_neighs(kk,2))==1
                                sum_neigh_targ_A=sum_neigh_targ_A+1;
                                
                            elseif A(targ_site_neighs(kk,1),targ_site_neighs(kk,2))==2
                                sum_neigh_targ_B=sum_neigh_targ_B+1;
                            end
                        end
                        
                        if neighbours_agent(kk,1) ~= 0 && neighbours_agent(kk,2) ~= 0 && neighbours_agent(kk,1) ~= nrows+1 && neighbours_agent(kk,2) ~= ncols+1
                            if A(neighbours_agent(kk,1),neighbours_agent(kk,2))==1
                                sum_neigh_agent_A=sum_neigh_agent_A+1;
                                
                            elseif A(neighbours_agent(kk,1),neighbours_agent(kk,2))==2
                                sum_neigh_agent_B=sum_neigh_agent_B+1;
                            end
                        end
                    end
                    
                    %movement probability of the focal agent
                    move_prob_agent=(1-p)^sum_neigh_agent_A*(1-r)^sum_neigh_agent_B;
                    
                    %calculate the movement probability of target agent
                    if A(target_site(1),target_site(2))==1
                        move_prob_targ=(1-p)^sum_neigh_targ_A*(1-r)^sum_neigh_targ_B;
                        
                    elseif A(target_site(1),target_site(2))==2
                        move_prob_targ=(1-q)^sum_neigh_targ_B*(1-r)^sum_neigh_targ_A;
                    end
                    
                    %generate two random numbers to check for the swap
                    r1=rand;
                    r2=rand;
                    
                    if r1<=move_prob_agent && r2<=move_prob_targ
                        %do swap
                        agent_val_copy=A(agent(1),agent(2));
                        targ_val_copy=A(target_site(1),target_site(2));
                        
                        A(target_site(1),target_site(2))=agent_val_copy;
                        A(agent(1),agent(2))=targ_val_copy;
                        
                    end
                    
                    %check for swap of type B agent
                elseif rnd<=rho && A(agent(1),agent(2))==2
                    
                    %find sum of neighbour of focal agent
                    sum_neigh_agent_A=0;
                    sum_neigh_agent_B=0;
                    
                    %find sum of neighbour of agent at target site
                    sum_neigh_targ_A=0;
                    sum_neigh_targ_B=0;
                    
                    %target site neighbours
                    targ_site_neighs=[target_site(1)-1 target_site(2);target_site(1)+1 target_site(2);...
                        target_site(1) target_site(2)-1;target_site(1) target_site(2)+1];
                    
                    %loop through each target site neighbour and
                    %calculate the number of neighbouring type A and
                    %type B agents
                    for kk=1:4
                        if targ_site_neighs(kk,1) ~= 0 && targ_site_neighs(kk,2) ~= 0 && targ_site_neighs(kk,1) ~= nrows+1 && targ_site_neighs(kk,2) ~= ncols+1
                            if A(targ_site_neighs(kk,1),targ_site_neighs(kk,2))==1
                                sum_neigh_targ_A=sum_neigh_targ_A+1;
                            elseif A(targ_site_neighs(kk,1),targ_site_neighs(kk,2))==2
                                sum_neigh_targ_B=sum_neigh_targ_B+1;
                            end
                        end
                        
                        if neighbours_agent(kk,1) ~= 0 && neighbours_agent(kk,2) ~= 0 && neighbours_agent(kk,1) ~= nrows+1 && neighbours_agent(kk,2) ~= ncols+1
                            if A(neighbours_agent(kk,1),neighbours_agent(kk,2))==1
                                sum_neigh_agent_A=sum_neigh_agent_A+1;
                                
                            elseif A(neighbours_agent(kk,1),neighbours_agent(kk,2))==2
                                sum_neigh_agent_B=sum_neigh_agent_B+1;
                            end
                        end
                    end
                    
                    %movement probability of the focal agent
                    move_prob_agent=(1-q)^sum_neigh_agent_B*(1-r)^sum_neigh_agent_A;
                    
                    %movement probability of the target site agent
                    if A(target_site(1),target_site(2))==1
                        move_prob_targ=(1-p)^sum_neigh_targ_A*(1-r)^sum_neigh_targ_B;
                    else
                        move_prob_targ=(1-q)^sum_neigh_targ_B*(1-r)^sum_neigh_targ_A;
                    end
                    
                    %generate two random number of check for swap
                    %between the focal agent and the target site agent
                    r1=rand;
                    r2=rand;
                    
                    if r1<=move_prob_agent && r2<=move_prob_targ
                        %do swap
                        agent_val_copy=A(agent(1),agent(2));
                        targ_val_copy=A(target_site(1),target_site(2));
                        
                        A(target_site(1),target_site(2))=agent_val_copy;
                        A(agent(1),agent(2))=targ_val_copy;
                        
                    end
                end
            end
                
        end
        
        %record the domain matix and video (if video==1)
        if t>=rec_times(j)
            rec_mat(:,:,j)=A;
            
            if video==1
                clf('reset');
                
                imagesc(A);

                colormap(c);
                axis square
                pbaspect([1 1 1]);
                title_ = sprintf('T = %f',t);
                title(title_);
                
                frame = getframe(gcf);
                writeVideo(v,frame);
                
            end
            
            %increment recording index
            j=j+1;
        end
        
    end
    
    rec_mat_full(:,:,:,i)=rec_mat;
    
end

%video
if video==1
    close(v);
end


% save workspace
if Save==1
    if IC==1
        save("two_dimensional_ts_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho)+".mat",'-v7.3')

    elseif IC==3 || IC==4

        save("pattern_abm_2D_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho)+".mat",'-v7.3')
    end

end