#include<bits/stdc++.h>

using namespace std;

const float RHO   = 0.20;
const float QNOT  = 0.99;
const float BETA  = 2.00;
const float ALPHA = 1.00;
const float GAMMA = 0.90;

const int MAX_ITERATIONS  = 100;
const int TOTAL_PHEROMONE = 30000;

const int RUN_COUNT  = 5;
const int TIME_LIMIT = 10;

// Ant
struct Ant{
  int currentGroup;
  set<int>visitedGroups;
  set<int>allowedPatients;
};

// Edge
struct Edge{
  bool  valid;
  float pheromone;
  float satisfactionRatio;
  Edge(){
    pheromone=0.0;
  }
};

// Group of patients
struct Group{
  int id;
  int patientCount;
  float satisfaction;
  set<int>patients;
};

//Node For random Selection
struct Node{
  int currentGroup;
  set<int>visitedGroups;
  set<int>allowedPatients;
};

// function checks mutual exclusivity in two groups
bool isMutuallyExclusive(set<int>&grp1,set<int>&grp2){
  for(auto i: grp1){
      if(grp2.find(i)!=grp2.end()){
        return false;
      }
    }
    return true;
}

// function checks wheather grp1 is a subset of grp2
bool isSubset(set<int>&grp1,set<int>&grp2){
  for(auto i: grp1){
    if(grp2.find(i)==grp2.end()){
      return false;
    }
  }
  return true;
}

// limits the total pheromone in matrix
void normalizePheromone(vector<vector<Edge>>&matrix){
  float totalpheromone=0.0;
  for(int i=0;i<matrix.size();i++){
    for(int j=0;j<matrix.size();j++){
      totalpheromone+= matrix[i][j].pheromone;
    }
  }
  for(int i=0;i<matrix.size();i++){
    for(int j=0;j<matrix.size();j++){
      matrix[i][j].pheromone *= (TOTAL_PHEROMONE/totalpheromone);
    }
  }
}

// calculating satis ratio for an edge
void fillSatisfactionRatio(vector<vector<Edge>>&matrix,vector<Group>&groups){
  for(int i=0;i<matrix.size();i++){
    float totalSatisfaction = 0.0;
    for(int j=0;j<matrix.size();j++){
      if(matrix[i][j].valid){
            totalSatisfaction+=groups[j].satisfaction;
      }
    }
    for(int j=0;j<matrix.size();j++){
      if(matrix[i][j].valid){
        matrix[i][j].satisfactionRatio=groups[j].satisfaction/totalSatisfaction;
      }
    }
  }
}

// finds next group randomly
int getNextRandomGroup(int currentGroup,set<int>&allowed,vector<vector<Edge>>&matrix,vector<Group>&groups){
  vector<int>options;
  for(int i=0;i<groups.size();i++){
    if(isSubset(groups[i].patients,allowed)){
      options.push_back(i);
    }
  }
  return options[rand()%options.size()];
}

// finds the next group for an ant
pair<int,float> getNextGroup(int currentGroup,set<int>&allowed,vector<vector<Edge>>&matrix,vector<Group>&groups){

  float q      = (rand()%100)/100.0;
  float maxx   = -100.0;
  float deltou = 0;
  int   nextGroup = -1;

  // if q <= qnot
  // then maximum of pheromone and satis ratio product is consider
  if(q<=QNOT){
    for(int i=0;i<groups.size();i++){
      if(isSubset(groups[i].patients,allowed) and pow(matrix[currentGroup][i].pheromone,ALPHA)*pow(matrix[currentGroup][i].satisfactionRatio,BETA) >maxx){
        maxx     = pow(matrix[currentGroup][i].pheromone,ALPHA)*pow(matrix[currentGroup][i].satisfactionRatio,BETA);
        nextGroup= i;
        deltou   = GAMMA*matrix[currentGroup][i].pheromone;
      }
    }
  }

  // if q > qnot then next group is selected from the probability function
  else{
    float prob=0.0,probSum=0.0;
    for(int i=0;i<groups.size();i++){
      if(isSubset(groups[i].patients,allowed)){
        probSum += pow(matrix[currentGroup][i].pheromone, ALPHA)*pow(matrix[currentGroup][i].satisfactionRatio, BETA);
      }
    }
    for(int i=0;i<groups.size();i++){
      if(isSubset(groups[i].patients,allowed)){
        prob += pow(matrix[currentGroup][i].pheromone, ALPHA)*pow(matrix[currentGroup][i].satisfactionRatio, BETA)/probSum;
        if(prob>=q){
          nextGroup=i;
          break;
        }
      }
    }
  }
  return {nextGroup,deltou};
}

// whenever an ant moves from one group to another the pheromone of the edge is updated
void localUpdate(vector<vector<Edge>>&matrix,int i,int j,float deltou){
  matrix[i][j].pheromone=(1-RHO)*matrix[i][j].pheromone + RHO*deltou;
  normalizePheromone(matrix);
}

// after an ant completes its iteration or finds an optimal solution the pheromone is updated on every edge
void globalUpdate(vector<Group>&groups,vector<vector<Edge>>&matrix,set<int>&solution){
  float deltapheromone;
  float maxpheromone = 0;

  //finding max pheromonefrom all the edges
  for(int i=0;i<groups.size();i++){
    for(int j=0;j<groups.size();j++){
      if(matrix[i][j].valid){
        if(matrix[i][j].pheromone>maxpheromone){
          maxpheromone = matrix[i][j].pheromone;
        }
      }
    }
  }
  deltapheromone = GAMMA*maxpheromone;
  for(int i=0;i<groups.size();i++){
    for(int j=0;j<groups.size();j++){

      //updating pheromone on each edge
      if(matrix[i][j].valid){
        matrix[i][j].pheromone = (1-RHO)*matrix[i][j].pheromone;
      }

      //  if edge is in optimaal solution it is increased little more by deltou
      if(solution.find(i)!=solution.end() && solution.find(j)!=solution.end()){
        if(matrix[i][j].valid){
          matrix[i][j].pheromone+=deltapheromone;
        }
      }
    }
  }
  normalizePheromone(matrix);
}


// including edges between mutually exclusive groups and initializing the pheromone
vector<vector<Edge>> constructGraph(vector<Group>&groups){
  vector<vector<Edge>>matrix(groups.size(),vector<Edge>(groups.size()));
  for(int i=1;i<groups.size();i++){
    for(int j=0;j<i;j++){
      if(isMutuallyExclusive(groups[i].patients,groups[j].patients)){
        matrix[i][j].valid     = matrix[j][i].valid     = true;
        matrix[i][j].pheromone = matrix[j][i].pheromone = 1.0;
      }
    }
  }

  // calling normaliztion of pheromone
  normalizePheromone(matrix);

  //filling the satisfaction ratios
  fillSatisfactionRatio(matrix,groups);

  return matrix;
}

// Random Selection Algorithm
set<int> randomSelection(Node &node,vector<vector<Edge>>&matrix,vector<Group>&groups){

  //remove patients of the group visited by ant
  for(int rem:groups[node.currentGroup].patients){
    node.allowedPatients.erase(rem);
  }
  node.visitedGroups.insert(node.currentGroup);

  // return groups set (solution when all the patients are covered
  if(node.allowedPatients.size()==0){
    return node.visitedGroups;
  }

  int nextGroup=getNextRandomGroup(node.currentGroup,node.allowedPatients,matrix,groups);

  if(nextGroup==-1){

    //if no group found
    node.visitedGroups.clear();
    return node.visitedGroups;
  }
  node.currentGroup=nextGroup;
  return randomSelection(node,matrix,groups);

}

// Ant Colony Optimization Algorithm
set<int> antColonyOptimization(Ant &ant,vector<vector<Edge>>&matrix,vector<Group>&groups){

  //remove patients of the group visited by ant
  for(int rem:groups[ant.currentGroup].patients){
    ant.allowedPatients.erase(rem);
  }
  ant.visitedGroups.insert(ant.currentGroup);

  // return groups set (solution when all the patients are covered
  if(ant.allowedPatients.size()==0){
    return ant.visitedGroups;
  }

  // calling select group function to find next group for ant to visit
  pair<int,float>p= getNextGroup(ant.currentGroup,ant.allowedPatients,matrix,groups);
  int   nextGroup = p.first;
  float delTou    = p.second;
  if(nextGroup==-1){

    // no group found return empty set
    ant.visitedGroups.clear();
    return ant.visitedGroups;
  }
  else{

    // update pheromone of the edges visited by ant
    localUpdate(matrix,ant.currentGroup,nextGroup,delTou);
    ant.currentGroup=nextGroup;

    // calling recursively aco function
    return antColonyOptimization(ant,matrix,groups);
  }
}

// local search optimizing the solution foound by the ant
set<int> localSearchRandomRestructure(Ant &ant,vector<vector<Edge>>&matrix,vector<Group>&groups,set<int>&solution){
  int currentGroup,r=rand()%solution.size();
  if(r==0)
    r++;

  //clearing certain part of the solution
  while(r--){
    currentGroup=*(solution.rbegin());
    solution.erase(currentGroup);
    ant.visitedGroups.erase(currentGroup);
    for(int j:groups[currentGroup].patients){
      ant.allowedPatients.insert(j);
    }
  }
  ant.currentGroup=currentGroup;

  //calling aco to reconstruct the solution
  return antColonyOptimization(ant,matrix,groups);
}

// find total satisfaction of the solution set of patients groups
float getSatisfaction(set<int>&solutionSet,vector<Group>&groups){
  float satisfaction=0.0;
  for(int i:solutionSet){
    satisfaction += groups[i].satisfaction;
  }
  return satisfaction;
}

float getUsingRandomSelection(vector<Group>&groups,set<int>&patients){

  clock_t start,end;

  start=clock();

  set<int> optimalSolution;
  float optimalSatisfaction = 0.0;

  vector<vector<Edge>>matrix=constructGraph(groups);

  // iterating over every ant
  while(true){
    set<int> solution;
    Node node;
    node.currentGroup    = rand()%groups.size();
    node.allowedPatients = patients;
    node.visitedGroups.clear();
    // calling aco
      solution = randomSelection(node,matrix,groups);

    // calculating satisfaction of the solution
    float satisfaction = getSatisfaction(solution,groups);

    // if the satisfaction is greater than optimal satisfaction so far
    if(satisfaction>optimalSatisfaction){
      optimalSolution = solution;
      optimalSatisfaction = satisfaction;
    }
    end=clock();
    if(double(end-start)/double(CLOCKS_PER_SEC)>=TIME_LIMIT){
    	return optimalSatisfaction;
    }
  }

  return optimalSatisfaction;
}

float getUsingAntColonyOptimization(vector<Group>&groups,set<int>&patients,bool localSearchFlag=false){

  clock_t start,end;

  start=clock();

  set<int> optimalSolution;
  float optimalSatisfaction = 0.0;

  vector<vector<Edge>>matrix=constructGraph(groups);

  // iterating over every ant
  while(true){
    set<int> solution;
    Ant ant;
    ant.currentGroup    = rand()%groups.size();
    ant.allowedPatients = patients;
    ant.visitedGroups.clear();

    int itCount = MAX_ITERATIONS;
    while(itCount--){

      // calling aco
      solution = antColonyOptimization(ant,matrix,groups);

      // calculating satisfaction of the solution
      float satisfaction = getSatisfaction(solution,groups);

      // if its zero then continuing
      if(satisfaction==0.0){
        continue;
      }

      // if the satisfaction is greater than optimal satisfaction so far
      if(satisfaction>optimalSatisfaction){
            optimalSolution = solution;
            optimalSatisfaction = satisfaction;
            globalUpdate(groups,matrix,optimalSolution);
      }

      // applying local search if flag is set
      if(localSearchFlag){
        solution         = localSearchRandomRestructure(ant,matrix,groups,solution);
        satisfaction     = getSatisfaction(solution,groups);
        if(satisfaction>optimalSatisfaction){
          optimalSolution = solution;
          optimalSatisfaction = satisfaction;
          globalUpdate(groups,matrix,optimalSolution);
        }
      }
      break;
    }
    end=clock();
    if(double(end-start)/double(CLOCKS_PER_SEC)>=TIME_LIMIT){
      return optimalSatisfaction;
    }
  }

  return optimalSatisfaction;
}


int main(){

  srand(time(0));

  int n;

  //input
  cin>>n;
  set<int> patients;
  vector<Group>groups(n);
  for(int i=0;i<n;i++){
    groups[i].id=i;
    cin>>groups[i].patientCount;
    cin>>groups[i].satisfaction;
    for(int j=0;j<groups[i].patientCount;j++){
      int patient;
      cin>>patient;
      patients.insert(patient);
      groups[i].patients.insert(patient);
    }
  }

  float totalSatisfaction=0.0;
  for(int run=1;run<=RUN_COUNT;run++){
    cout<<endl<<"Run "<<run<<" Using Random Selection ";
    float satisfaction = getUsingRandomSelection(groups,patients);
    cout<<endl<<"Maximum Satisfaction Found : "<<satisfaction<<endl;
    totalSatisfaction+=satisfaction;
  }
  cout<<endl<<"Average satisfaction : "<<totalSatisfaction/RUN_COUNT<<endl;

  cout<<endl<<"----------------------------------------------------------------------------------"<<endl;

  totalSatisfaction=0.0;
  for(int run=1;run<=RUN_COUNT;run++){
    cout<<endl<<"Run "<<run<<" Using ACO Without LS ";
    float satisfaction = getUsingAntColonyOptimization(groups,patients);
    cout<<endl<<"Maximum Satisfaction Found : "<<satisfaction<<endl;
    totalSatisfaction+=satisfaction;
  }
  cout<<endl<<"Average satisfaction : "<<totalSatisfaction/RUN_COUNT<<endl;


  cout<<endl<<"----------------------------------------------------------------------------------"<<endl;

  totalSatisfaction=0.0;
  for(int run=1;run<=RUN_COUNT;run++){
    cout<<endl<<"Run "<<run<<" Using ACO With Random Restructure Local Search ";
    float satisfaction = getUsingAntColonyOptimization(groups,patients,true);
    cout<<endl<<"Maximum Satisfaction Found : "<<satisfaction<<endl;
    totalSatisfaction+=satisfaction;
  }
  cout<<endl<<"Average satisfaction : "<<totalSatisfaction/RUN_COUNT<<endl;

  cout<<endl<<"----------------------------------------------------------------------------------"<<endl;

  cout<<endl;
  return 0;
}
