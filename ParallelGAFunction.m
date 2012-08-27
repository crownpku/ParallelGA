function ParallelGAFunction()
clear
%%Parameters for Wiestrass Problem
%KKK=20;%Number of k
%DDD=1;%Number of Xi
%NumberArea = 15; %The number of locus to represent one number



%%Parameters for one population
nLocus = 14;
nChromosome = 10;
nNode = 1;
f = zeros(nNode, nChromosome); %%The fitness of the chromosomes
o = zeros(nNode, nChromosome); %%The fitness order of the chromosomes in one node
pr = zeros(nNode, nChromosome); %%Row Probility to survive
pc = zeros(nNode, nLocus); %%Locus Probility to mutate
N = zeros(nLocus, 2);%%Count, 1means0, 2means1
Np = zeros(nLocus, 2);%%Weighter Count, 1means0, 2means1

%%Parameters for whole problem
AdjMatrix = zeros (nNode, nNode);
Population = zeros (nNode, nChromosome, nLocus);
CommunicationPeriod = 1; %Number of Generations before communication
RoundNum = 1000; %Number of total generations = CommunicationPeriod * RoundNum
CommuChromosomeNum=1;%Number of exchange chromosome when communicating
runtimes=100; %For each cap, run runtimes and take average




%%Function to initialize
    %{
      function init(Node)
      end
%}


%%Function to Calculate the fitness
    function [Fitness] = fitness( Chromosome )
        x=0;
        y=0;
        for i=1:(nLocus/2)
            x = x+ (2^i) * Chromosome(i);
            y = y+ (2^i) * Chromosome(i+ nLocus/2);
        end
        x = -5 + x*10/128;
        y = -5 + y*10/128;
        %Fitness = -100*x*sin(abs(100*x)^(0.5))-100*y*sin(abs(100*y)^(0.5));
        %Fitness = x^2+y^2-10*(cos(2*pi*x)+cos(2*pi*y))+20;
        Fitness = -20*exp(-0.2* sqrt( (x*x+y*y)/2) ) - exp( (cos(2*pi*x)+cos(2*pi*y))/2    )+20+exp(1);
    end


%%Function to sort the fitness
    function sortPop (Node, a, b)
        i=a; j=b;
        m=f(Node, o(Node, ceil((a+b)/2)) );
        while i<=j
            while f(Node, o(Node, i) )<m
                i=i+1;
            end
            while f(Node, o(Node, j) )> m
                j=j-1;
            end
            if i<=j
                tmp=o(Node, i);
                o(Node, i)=o(Node, j);
                o(Node,j) = tmp;
                i=i+1; j=j-1;
            end
        end
        if i<b
            sortPop(Node, i,b);
        end
        if j>a
            sortPop(Node, a,j);
        end
    end

%%Function to generate a new chromosome
    function [Chromosome] = genChromosome (Node, Chromosome )
        for i = 1:nLocus
            if rand()<pc(Node, i)
                Chromosome(i)=1-Chromosome(i);
            end
        end
    end

%%Function to generate the next generation
    function genNextGeneration(Node)
        for i=1:nChromosome
            if rand()>pr(Node, i)
                %Matlab only allows value transfer
                Population(Node, o(Node, i),: )=genChromosome(Node, Population(Node, o(Node, i),: ) );
                f(Node, o(Node, i))= fitness (Population(Node, o(Node, i),: ) );
            end
        end
        for i = 1:nChromosome
            o(Node, i)=i;
        end
        sortPop(Node, 1, nChromosome);
    end

%%Fcuntion to make the first generation
    function firstGeneration(Node)
        for i=1:nChromosome
            o(Node,i)=i;
        end
        Population(Node,:,:) = zeros (nChromosome, nLocus);
        for i=1:nChromosome
            pr(Node,i)=0;
        end
        for i=1:nLocus
            pc(Node, i)=0.5;
        end
        genNextGeneration(Node);
        for i=1:nChromosome
            pr(Node, i)=1- (i-1)/nChromosome;
        end
    end
%%Fucntion to coleect statistics and calculate pc
    function collectStat(Node)
        N = zeros(nLocus, 2);
        Np = zeros(nLocus, 2);
        for i=1:nLocus
            for j=1:nChromosome
                N(i, Population(Node, o(Node, j),i )+1 )= N(i,Population(Node, o(Node, j),i )+1 )+1;
                Np(i, Population(Node, o(Node, j),i )+1)=Np(i, Population(Node, o(Node, j),i )+1)+pr(Node,j);
            end
        end
        for i=1:nLocus
            tmp1=0; tmp2=0;
            for j=1:2
                tmp1 = tmp1 + N(i, j);
                tmp2 = tmp2 + Np(i, j);
            end
            for j=1:2
                N(i ,j)=N(i ,j)/tmp1;
                Np(i,j)=Np(i, j)/tmp2;
            end
        end
        for i=1:nLocus
            pc(Node, i)=1- abs(Np(i,1)-0.5)*2;
        end
    end

%Funtion to make communications for all the populations
    function communicate()
        for j=1:nNode
            immicount=-1;
            for i=1:nNode
                if AdjMatrix(i,j)==1
                    immicount=immicount+1;%Node i Immigrate to Node j
                    for k=1:CommuChromosomeNum
                        Population(j, o(j, nChromosome-(k+immicount*CommuChromosomeNum)+1),:)=Population(i, o(i,k), :);
                        f(j,o(j,nChromosome-(k+immicount*CommuChromosomeNum)+1) ) = f(i, o(i,k));
                    end
                end
            end
        end
        
        
        %{
        for i=1:CommuChromosomeNum
            for j=1:nLocus
                Population(Node1, o(Node1,nChromosome-i+1), j)=Population(Node2, o(Node2,i), j);
                Population(Node2, o(Node2,nChromosome-i+1), j)=Population(Node1, o(Node1,i), j);
            end
            f(Node1, o(Node1,nChromosome-i+1))=f(Node2, o(Node2,i) );
            f(Node2, o(Node2,nChromosome-i+1))=f(Node1, o(Node1,i) );
        end
        for i = 1:nChromosome
            o(Node1, i)=i;
            o(Node2, i)=i;
        end
        sortPop(Node1, 1, nChromosome);
        sortPop(Node2, 1, nChromosome);
        %}
        
        
    end

%%The main function starts here.
%%First initialize the AdjMatrix.
    function main()
        %No Links
        %Pairs
        %{
        AdjMatrix(1,2)=1;
        AdjMatrix(2,1)=1;
        AdjMatrix(3,4)=1;
        AdjMatrix(4,3)=1;
        AdjMatrix(5,6)=1;
        AdjMatrix(6,5)=1;
        AdjMatrix(7,8)=1;
        AdjMatrix(8,7)=1;
       %}
        
        %Ring
        %{
        AdjMatrix(1,2)=1;AdjMatrix(1,8)=1;
        AdjMatrix(2,1)=1;AdjMatrix(2,3)=1;
        AdjMatrix(3,2)=1;AdjMatrix(3,4)=1;
        AdjMatrix(4,3)=1;AdjMatrix(4,5)=1;
        AdjMatrix(5,4)=1;AdjMatrix(5,6)=1;
        AdjMatrix(6,5)=1;AdjMatrix(6,7)=1;
        AdjMatrix(7,6)=1;AdjMatrix(7,8)=1;
        AdjMatrix(8,7)=1;AdjMatrix(8,1)=1;
        %}
        
        
        
        %XXXX Band
        %{
        AdjMatrix(1,2)=1;AdjMatrix(1,4)=1;AdjMatrix(1,6)=1;AdjMatrix(1,8)=1;
        AdjMatrix(2,1)=1;AdjMatrix(2,3)=1;AdjMatrix(2,5)=1;AdjMatrix(2,7)=1;
        AdjMatrix(3,2)=1;AdjMatrix(3,4)=1;AdjMatrix(3,6)=1;AdjMatrix(3,8)=1;
        AdjMatrix(4,3)=1;AdjMatrix(4,1)=1;AdjMatrix(4,7)=1;AdjMatrix(4,5)=1;
        AdjMatrix(5,6)=1;AdjMatrix(5,8)=1;AdjMatrix(5,2)=1;AdjMatrix(5,4)=1;
        AdjMatrix(6,5)=1;AdjMatrix(6,7)=1;AdjMatrix(6,1)=1;AdjMatrix(6,3)=1;
        AdjMatrix(7,6)=1;AdjMatrix(7,8)=1;AdjMatrix(7,2)=1;AdjMatrix(7,4)=1;
        AdjMatrix(8,7)=1;AdjMatrix(8,5)=1;AdjMatrix(8,3)=1;AdjMatrix(8,1)=1;
        %}
        
        %Fully Connected
        %{
        for i=1:nNode
            for j=1:nNode
                if i~=j
                    AdjMatrix(i,j)=1;
                end
            end
        end
        %}
        
        
        
        fil = fopen ('Ackley10Chromosomes1000Generations100Times.txt', 'wt');
        
        %{
        for i=1:nNode
            init(i);
        end
        %}
        
        
        %Run runtimes and get the average result for each generation
        GeneAve(1:RoundNum) = 0;
        for roundcountt=1:runtimes
            for i=1:nNode
                firstGeneration(i);
                %Subspaces
                %{
                for j=1:nChromosome
                    if i==1
                        Population(i, j, 1)=0; Population(i, j, 2)=0; Population(i, j, 3)=0;
                    end
                    if i==2
                        Population(i, j, 1)=0; Population(i, j, 2)=0; Population(i, j, 3)=1;
                    end
                    if i==3
                        Population(i, j, 1)=0; Population(i, j, 2)=1; Population(i, j, 3)=0;
                    end
                    if i==4
                        Population(i, j, 1)=0; Population(i, j, 2)=1; Population(i, j, 3)=1;
                    end
                    if i==5
                        Population(i, j, 1)=1; Population(i, j, 2)=0; Population(i, j, 3)=0;
                    end
                    if i==6
                        Population(i, j, 1)=1; Population(i, j, 2)=0; Population(i, j, 3)=1;
                    end
                    if i==7
                        Population(i, j, 1)=1; Population(i, j, 2)=1; Population(i, j, 3)=0;
                    end
                    if i==8
                        Population(i, j, 1)=1; Population(i, j, 2)=1; Population(i, j, 3)=1;
                    end
                end
                %}
            end
            
            
            for i=1:RoundNum
                for j=1:CommunicationPeriod
                    for k=1:nNode
                        collectStat(k);
                        genNextGeneration(k);
                    end
                end
                communicate();
                for j=1:nNode
                    for k= 1:nChromosome
                        o(j, k)=k;
                    end
                    sortPop(j, 1, nChromosome);
                end
                 %%Find the best fitness from all the nodes
                BestFitness=1e10;
                for j=1:nNode
                    if f(j,o(j,1))<BestFitness
                        BestFitness = f(j,o(j,1));
                    end
                end
                GeneAve(i) = GeneAve(i) + BestFitness;
            end
        end
        %Calculate the Standard Deviation
        for roundcountt = 1:RoundNum
            fprintf(fil,'%d\t', roundcountt );
            fprintf(fil,'%f\t', GeneAve(roundcountt)/runtimes );
            fprintf(fil,'\n');
        end
        fclose(fil);
    end

main();


end