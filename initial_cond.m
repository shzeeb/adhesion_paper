%initial condition for the lattice model:
%nrows: number of rows
%ncols: number of columns
%type 1: Figures 4-6, 10 -- all sites in the middle 20% occupied by speices A. The rest
% occupied by speices B with density 0.5.
%type 2: Figure 1 -- single species - middle 20% occupied
%type 3: Figure 7 (d) and (g) uniformly at random: Species A initilaised at density 0.45,
%species B also at density of 0.45. Remaining 10% unccupied.
%type 4: Figure 7 (a) -- left half populated uniformly at random with species A at density
%0.9 and the right left with species B at density 0.9.
%type 5: Figure 11 (a) -- species A at density 0.7 uniformly at random and
%species B at density 0.3.
%type 6: Figure 11 (d) -- species A and B at densoty 0.5 uniformly at
%random.

function M=initial_cond(nrows,ncols,type)

if type==1
    %left matrix (50% initialised with species B)
    A=rand([nrows,floor(0.4*ncols)]);

    %middle matrix (all species A)
    B=ones(nrows,ceil(0.6*ncols)-floor(0.4*ncols));

    %right matrix (50% initialised with species B)
    C=rand([nrows,floor(0.4*ncols)]);
    
    A=2*(A<0.50);
    
    C=2*(C<0.50);
    

    %concaterate A, B and C
    M=[A,B,C];

    
elseif type==2
    A=zeros(nrows,ncols);
    
    A(:,ceil(0.45*ncols)+1:floor(0.55*ncols))=1;
    
    M=A;


elseif type==3

    M=rand(nrows,ncols);

    M(M<=0.45)=1;
    M(M>0.45 & M<=0.90)=2;
    M(M>0.90 & M<1)=0;


elseif type==4
    
    A=rand(nrows,ceil(ncols/2));
    B=rand(nrows,ceil(ncols/2));


    A(A<=0.90)=1;
    A(A>0.90 & A<1)=0;
    B(B<=0.90)=2;
    B(B>0.90 & B<1)=0;

    M=[A B];

elseif type==5

    M=rand(nrows,ncols);

    M(M<=0.70)=1;
    M(M>0.70)=2;

elseif type==6

    M=rand(nrows,ncols);

    M(M<=0.50)=1;
    M(M>0.50)=2;
    
end

end