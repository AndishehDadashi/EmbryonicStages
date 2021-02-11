#' Author and Code Owner: Andisheh Dadashi, @AndishehDadashi
#'
#' Date: April 16th, 2020
#'
#' Title: Network Expansion
#'
#' Description: Constructing and expanding the metabolic network of
#' the embryonic stages starting from an initial seed (set of essential metabolites).
#'
#' Human 1 Model used for this study can be accessible from https://github.com/SysBioChalmers/Human-GEM
#'
#' @details:
#'
#' In this program we defined three functions and one method as follow:
#' 1. Function LP: LP function is a linear programming function that requires the installation of
#'    Gurobi solver on the system. (https://www.gurobi.com/downloads/)
#'    LP requires the information regarding the linear function and constraints.
#' 2. Function Adjusting_Uptake: Using this function we are able to adjust the secretion
#'    and uptake rate of the metabolites base on the stage.
#' 3. Function Reconstruct: Reconstruct function changes the reaction set for each stage
#'    and creates the seed set,  the matrix of products P and reactants R for each stage
#' 4. Method Expansion: Expanding the stages

#' ########################################
#' ###                                  ###
#' ### To load the required libraries   ###
#' ###                                  ###
#' ########################################

library("data.table", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("ggplot2", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("stringr", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("lpSolve", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("ggpubr", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("Matrix", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("dplyr", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("gdata", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("slam", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("tools",lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("pryr", lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("grid",lib.loc='~/R/x86_64-pc-linux-gnu-library/3.6');
library("gurobi");

# Helper function for concatenating the strings

`+` <- function(e1, e2) {
	if (is.character(e1) | is.character(e2)) {
		paste0(e1, e2)
	} else {
		base::`+`(e1, e2)
	}
}

#' ########################################
#' ###                                  ###
#' ###        Reading the inputs        ###
#' ###                                  ###
#' ########################################
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!= 2) {
	stop("Usage: <Output_ES Result Path> <Data Directory Path>", call.=FALSE)
}

# Assign the paths from arguments
Output_ES<-args[1]
Data<-args[2]

# Create the Output directory if it does not already exist
if (!file.exists(Output_ES)){
	dir.create(Output_ES)
}

Environment<-read.csv(Data + "Environment.csv",sep =",",header = TRUE);
smat<-read.table(Data + "Smat.csv" ,sep =",", header = FALSE);
Metabolites_org<-read.csv(Data + "Metabolites.csv",sep =",",header = TRUE); # Original Metabolite file
Reaction_Info_File<-read.csv(Data + "Reactions.csv",sep =",",header = TRUE); # read the reaction information file
Stage_Reactions_file<-read.csv(Data + "AllStages.csv",sep =",",header = FALSE); # contains the necessary reactions for each stage
assign("Stage_Reactions_file",Stage_Reactions_file,envir = .GlobalEnv )
NoENSGReactions_file<-read.csv( Data + "NoENSGReactions.csv",sep =",",header = FALSE); # the reaction list without ENSG index

# smat is the universe stoichiometric matrix. Rows are Metabolites and columns are reactions
smat<-as.matrix(smat);
# Column of abbreviation of metabolites
num_cols<-ncol(smat)
num_rows<-nrow(smat)
Reactions_org<-as.character(Reaction_Info_File[,1]); # Change the type to character
Metabolites_org<-as.character(Metabolites_org$Abb) # Change the type to character
All_obj<-c() # making vectors to save the objectives for each generated network
Rxnsize_vect<-c() # making vectors to save the number of reactions before the non-viability
metsize_vect<-c() # making vectors to save the number of metabolites before the non-viability
Generations<-c()
# Random_Metabolite_name<-c() # making vectors to save the random metabolites added to each network
# matrix to save information in the generated network
assign("vect_1", matrix(rep(" ",(num_cols+1)*2),nrow=(num_cols+1), ncol=2 ), envir = .GlobalEnv )
assign("Metabolites_Generation",matrix(rep(" ",num_rows*num_rows),nrow=num_rows, ncol=num_rows)  , envir = .GlobalEnv)
assign("Reactions_Flux_Generation",matrix(rep(" ",num_rows*num_rows),nrow=num_rows, ncol=num_rows) , envir = .GlobalEnv)
assign("Reactions_Flux_RealFormula",matrix(rep(" ",num_rows*num_rows),nrow=num_rows, ncol=num_rows) , envir = .GlobalEnv)

Total_stages<- 5 # Counting stages beginning from the first stage
Generation <- 1 # set the initial generation to one
Core_Reactions<-c() # changing the name of reactions set
Core_Metabolites<-c() # changing the name of seed metabolite
List_1<-list()
List_1[[1]]<- c("Stage","sub_Generation", "Number of generation","Growth values", "Number of Metabolites" , "Number of Reactions" )
count_Info<-3 # Counter
count_reaction_flux<-1 # Counter


#' ################################
#' ###                          ###
#' ###    Solving LP For FBA    ###
#' ###                          ###
#' ################################

LP<-function(Smat, Objective,rhs, rhseqs, lb, ub){
	# program maximizes
	model<-list();
	model$modelsense<-"max";
	model$vtype<-"C";
	model$obj<-Objective;
	model$A<-Smat;
	model$rhs<-rhs;
	model$sense<-rhseqs;
	model$lb<-lb;
	model$ub<-ub;
	params<-list(Method=1,OutputFlag=1,OptimalityTol=1e-4)
	# Method 3 for LP (dual simplex) and Method 1 when memory is limited
	solution<-gurobi(model, params=params); # program maximizes
	rm(model, params)
	return(solution)
}

#' ################################
#' ###                          ###
#' ###    Adjusting_Uptake      ###
#' ###                          ###
#' ################################

Adjusting_Uptake<-function(stage_counter,Reaction_list,lb){
	# we can set up uptake and secretion rate of any metabolites
	Pyruvate_Carbon <- which( Reaction_list == "HMR_9133");
	Glucose_Carbon <- which( Reaction_list == "HMR_9034");
	# To change the uptake (lb) and secretion (ub) rate in different stages
	# if(stage_counter <= 2){
	#   lb[Glucose_Carbon]<- -2
	#   lb[Pyruvate_Carbon]<- -1
	# }else if(stage_counter <= 4){ # 2 cell and 4 cell
	#   lb[Glucose_Carbon]<- -5
	#   lb[Pyruvate_Carbon]<- -3
	# }else if(stage_counter == 5){ # 8 cell
	#   lb[Glucose_Carbon]<- -2.59
	#   lb[Pyruvate_Carbon]<- -5.43
	# }
	lb<-as.numeric(lb); # to make sure this is going to be numeric values
	assign("lb",lb, envir = .GlobalEnv) # save to the GlobalEnv
}

#' #############################################
#' ###                                       ###
#' ###  Reconstruct the seeds for each stage ###
#' ###                                       ###
#' #############################################

Reconstruct<-function(smat,stage_counter,Reactions,Reactions_id,Metabolites,Compounds_id,Reaction_Info_File,NoENSGReactions_file){
	
	Stage_Reactions<-Stage_Reactions_file[,stage_counter] # Reading the Stage_Reactions file
	Stage_Reactions_indices<-which(Reactions %in% Stage_Reactions); # Finding the index of reactions in the Stage_Reactions
	NoENSGReactions<-NoENSGReactions_file[,1]
	NoENSGReactions_indices<-which(Reactions %in% NoENSGReactions);
	# from previous stage keep the ones that are expressed in the new stage
	Gene_difference<-Reactions_id[which(Reactions_id %in% Stage_Reactions_indices)]
	ATPS4m<-which( Reactions_org  == "HMR_6916") # ATPS4m reaction
	# Recognizing and organizing the objective functions
	Biomass_Human<-which( Reactions_org  == "biomass_human")
	Biomass_biomass_Recon3D<-which( Reactions_org  == "biomass_Recon3D")
	Biomass_maintenance_Recon3D<-which( Reactions_org == "biomass_maintenance_Recon3D")
	Biomass_maintenance_noTrTr_Recon3D<- which( Reactions_org  == "biomass_maintenance_noTrTr_Recon3D")
	All_Biomass<-c(Biomass_Human,Biomass_biomass_Recon3D,Biomass_maintenance_Recon3D,
	Biomass_maintenance_Recon3D,Biomass_maintenance_noTrTr_Recon3D)
	ToBeKept_Reactions<-unique(c(All_Biomass,ATPS4m, Gene_difference,
	Stage_Reactions_indices, NoENSGReactions_indices))
	assign("ToBeKept_Reactions",ToBeKept_Reactions, envir = .GlobalEnv) # save to the GlobalEnv
	Reconstructed_Smat<-smat[,ToBeKept_Reactions] # Keeping the reaction columns we need
	num_rows<-nrow(Reconstructed_Smat); # number of rows
	num_cols<-ncol(Reconstructed_Smat); # number of columns
	d<-num_rows*num_cols  # dimension of smat
	# Make the matrix of products P and the matrix of reactants R
	P<-matrix(rep(0,d),nrow = num_rows,ncol =num_cols  )
	R<-matrix(rep(0,d),nrow = num_rows,ncol =num_cols )
	for(i in 1:num_rows){
		for(j in 1:num_cols){
			if( Reconstructed_Smat[i,j] > 0 ){P[i,j]<- 1}
			if( Reconstructed_Smat[i,j] < 0 ){R[i,j]<- 1}
		}
	}
	# Make sparse matrix
	P<-Matrix(P[,], sparse=TRUE);
	R<-Matrix(R[,], sparse=TRUE);
	# Add all the consumed metabolites in each reaction and count them.
	B<-matrix(rep(0,num_cols),nrow = num_cols ,ncol=1 )
	for(j in 1:num_cols) {
		B[j,1]<-sum(R[,j])
	}
	Reconstructed_ReactionList<-Reaction_Info_File[ToBeKept_Reactions,]
	lb<-Reconstructed_ReactionList$LB; # updating the uptake values of the reactions
	ub<-Reconstructed_ReactionList$UB; # updating the secretion values of the reactions
	Objective<-Reconstructed_ReactionList$c; # Update the objective
	# another copy of Reaction file so I can modify it and still have the original one
	# column of all reactions in reaction file
	Reconstructed_Reaction<-as.character(Reconstructed_ReactionList[,1]); # Change the type to character
	# change their ids base on the new list
	Reconstructed_Reactions_id<-unique(which(Reconstructed_Reaction %in% Reactions[Reactions_id]))
	# save them to the GlobalEnv
	assign("P",P, envir = .GlobalEnv)
	assign("R",R, envir = .GlobalEnv)
	assign("B",B, envir = .GlobalEnv)
	assign("Reconstructed_Smat",Reconstructed_Smat, envir = .GlobalEnv)
	assign("Reconstructed_Reaction",Reconstructed_Reaction, envir = .GlobalEnv)
	assign("Objective",Objective, envir = .GlobalEnv)
	assign("ub",as.numeric(ub), envir = .GlobalEnv)
	assign("lb",as.numeric(lb), envir = .GlobalEnv)
	assign("Compounds_id",Compounds_id, envir = .GlobalEnv)
	assign("Reconstructed_Smat",Reconstructed_Smat, envir = .GlobalEnv)
	assign("Reconstructed_Reactions_id",Reconstructed_Reactions_id,envir = .GlobalEnv)
	print(c("Number of metabolites in the seed to begin with: ",length(Compounds_id) ),sep="\n",quote=FALSE);
	print(c("Number of reactions that are not silent: " ,length(ToBeKept_Reactions)),sep="\n",quote=FALSE);
}

# creating a function that reads all the inputs and generates the networks and calculates the growth rates

Create_Output<-function(filename,argument1) {
	Df<-data.frame(argument1)
	write.table(Df,Output_ES + filename ,sep = "," ,row.names = FALSE,
	col.names = FALSE, quote = FALSE);
}

Create_Output_cbind<-function(filename,argument1,argument2,argument3) {
	Df<-cbind(argument1,argument2,argument3)
	write.table(Df,Output_ES + filename ,sep = "," ,row.names = FALSE,
	col.names = FALSE, quote = FALSE);
}

#' ################################
#' ###                          ###
#' ### NetWork Expansion Method ###
#' ###                          ###
#' ################################

# stage_counter<-1
for( stage_counter in 1:Total_stages){
	
	print(c("Begining of Stage : ",stage_counter),sep="\n",quote=FALSE);
	# For each new environment reset the information
	Smat<-smat # Back to the original data for the new environment
	Reactions <- Reactions_org # making another copy of Reaction file
	Metabolites <- Metabolites_org # making another copy of metabolites
	Reactions_id<- 1 # To store the reaction's id numbers starting with biomass. The first column
	Mt<-as.character(Environment[,1]) # metabolites in the environment file
	Compounds_id<-c() # create an empty vector to store the metabolites id numbers in the main metabolite file
	for (i in 1:length(Mt)) {
		Compounds_id[i]<- which( Metabolites == Mt[i])
		# compare each metabolite in the environment with the metabolite file and find the id number for that metabolite
	}
	# Updating the seed set
	if( length(Core_Metabolites) > 0){Compounds_id<-unique(c(Compounds_id,Core_Metabolites))}
	if( length(Core_Reactions) > 0){Reactions_id <-unique(c(Reactions_id,Core_Reactions))}
	# lb[Reactions == "HMR_9048"]<- 0 # Anaerobic: set the Environment to block Oxygen
	
	# Using the function Reconstruct to create P R matrices for each stage
	Reconstruct(smat,stage_counter,Reactions, Reactions_id, Metabolites,Compounds_id, Reaction_Info_File,NoENSGReactions_file)
	
	num_rows<-nrow(Reconstructed_Smat); # number of rows
	num_cols<-ncol(Reconstructed_Smat); # number of columns
	rhs<- rep(0,num_rows); # Modify right hand side of equations
	rhseqs<-rep("=",num_rows);  # Modify the equality sign of of equations
	
	metsize <- length(Compounds_id) # number of metabolites in the seed
	metsize_int <- 0 # assign zero to the number of metabolite
	last_obj<-0;
	sub_Generation<-1
	solution<-list()
	
	# continue adding metabolites and reactions as long as there is a possibility to add more metabolites
	#   and the growth rate of objective function is zero
	while ( metsize > metsize_int ) {
		List_Metabolites<-c()
		List_Reactions<-c()
		List_Flux<-c()
		print(c("Generation: ",Generation),sep="\n",quote=FALSE);
		
		metsize_int<- metsize # set the number of metabolites to the seed set (update the previous seed set)
		X<-matrix(rep(0,num_rows*1),nrow =num_rows,ncol=1 ) # creating a matrix to store the added metabolites
		X[Compounds_id,1]<- 1 # assign one to those metabolites in the seed set
		Sum_Of_Coef<- t(R) %*% X
		# metabolites presented (sum of all the coefficient of metabolites in each reaction) in each reaction.
		Y<-matrix(rep(0,num_cols*1),nrow = num_cols,ncol=1 )
		# if Sum_Of_Coef is equal to the total number of consumed metabolites (reactants) in a reaction then choose
		#   that reaction to be added to the network
		for (j in 1:num_cols){if (Sum_Of_Coef[j] == B[j] ) {Y[j,1]<-1} else {Y[j,1]<-0}}
		All_Biomass<-1:4
		Core_Reactions<-unique(c(All_Biomass,which(Y == 1))) # add the reactions to the previous list
		
		X_Updated<-matrix(rep(0,num_rows*1),nrow =num_rows,ncol=1 )
		P_Y<- P %*% Y # Find the reaction which produces the metabolites.
		# if P_Y is positive for a metabolite it means that metabolite can be produced
		for (i in 1:num_rows){if (P_Y[i,1] >0 ) {X_Updated[i,1] <-1} else {X_Updated[i,1]<-0}}
		Core_Metabolites<-unique(c(Compounds_id,which(X_Updated == 1) ));  # add the metabolite to the previous list
		Compounds_id<- Core_Metabolites
		metsize <-length(Core_Metabolites) # number of metabolites in seed set
		# save to the GlobalEnv
		assign("Core_Reactions",Core_Reactions , envir = .GlobalEnv)
		assign("Core_Metabolites",Core_Metabolites, envir = .GlobalEnv)
		# create a new smat with these metabolites and reactions
		assign("Smat_Core",Reconstructed_Smat[Core_Metabolites, Core_Reactions] , envir = .GlobalEnv)
		
		# solve the new model (based on the new expanded network)
		#   program maximizes
		rhs_FBACheck<- rep(0,nrow(Smat_Core));
		rhseqs_FBACheck<-rep("=",nrow(Smat_Core) );
		Objective_FBACheck<-Objective[Core_Reactions];
		ub_FBACheck<-as.numeric(ub[Core_Reactions]);
		Adjusting_Uptake(stage_counter,Reconstructed_Reaction,lb);
		lb_FBACheck<-as.numeric(lb[Core_Reactions]);
		solution<-LP(Smat_Core, Objective_FBACheck ,rhs_FBACheck, rhseqs_FBACheck, lb_FBACheck, ub_FBACheck);
		assign("solution",solution, envir = .GlobalEnv);
		S <- solution[["status"]];
		# check for optimal solution. If not optimal set everything to zero to repeat the process again
		if(S == "OPTIMAL"){
			last_obj<-unlist(solution[6]); last_obj<-last_obj[[1]];
			# extracting flux of objective value from the solutions
			# if(abs(last_obj) >=  MaxFlux) { last_obj<- 0} # reset the objective value if it is not optimal
			last_obj<- assign("last_obj",last_obj, envir = .GlobalEnv);
			List_1[[count_Info]]<- c(stage_counter , sub_Generation, Generation,last_obj,length(Core_Metabolites),length(Core_Reactions) )
			count_Info<-count_Info+1;
			assign("List_1",List_1,envir = .GlobalEnv)
			Generation<-Generation+1 # update the Generation the initial point is one
			sub_Generation<-sub_Generation+1
			Generation<-assign("Generation",Generation, envir = .GlobalEnv)
		} else {last_obj<- assign("last_obj",0, envir = .GlobalEnv);  }
		# save to the GlobalEnv
		assign("sol_x",unlist(solution[7]), envir = .GlobalEnv) # save to the Global Env
		assign("last_obj",last_obj, envir = .GlobalEnv);  # save to the Global Env
		
		List_Metabolites<-c("Metabolite_" + stage_counter , "Generation_" +(Generation-1),Metabolites[Core_Metabolites])
		List_Reactions<-c("Reaction_" + stage_counter , "Generation_" +(Generation-1),Reconstructed_Reaction[Core_Reactions])
		List_Flux<-c("Flux_" + stage_counter , "Generation_" +(Generation-1) ,unlist(solution[7]) )
		Metabolites_Generation[1:length(List_Metabolites), (Generation-1)]<-List_Metabolites
		Reactions_Flux_Generation[1:length(List_Reactions),count_reaction_flux ]<-List_Reactions
		Reactions_Flux_Generation[1:length(List_Flux),count_reaction_flux+1]<-List_Flux
		Reactions_Flux_RealFormula[1:length(List_Reactions),count_reaction_flux ]<- c("Reaction_" + stage_counter , "Generation_" +(Generation-1),
		Reaction_Info_File[which(Reactions_org %in% List_Reactions),5])
		Reactions_Flux_RealFormula[1:length(List_Flux),count_reaction_flux+1]<-List_Flux
		count_reaction_flux<-count_reaction_flux+2
		
		# Create_Output("SubnetSmat" + (Generation-1)+ ".csv",Smat_Core)
		Create_Output_cbind("SubnetReactions_Info" + (Generation-1) + ".csv",Reconstructed_Reaction[Core_Reactions],lb_FBACheck,ub_FBACheck )
		Create_Output("SubnetMetabolites_Info" + (Generation-1)+ ".csv",Metabolites[Core_Metabolites])
		
	}
	
	assign("Sub_Generation_" + stage_counter, sub_Generation-1 , envir = .GlobalEnv) # save the Generation's number to Global Env
	assign("sol_obj_" + stage_counter ,last_obj, envir = .GlobalEnv) # save to the Global Env
	assign("Whole_solution_" + stage_counter , solution, envir = .GlobalEnv) # save the last solution to Global Env
	assign("Rxnsize_" + stage_counter,ncol(Smat_Core), envir = .GlobalEnv) # number of reactions in new s mat
	assign("Metsize_" + stage_counter,nrow(Smat_Core), envir = .GlobalEnv) # number of metabolites in new s mat
	assign("Sub_Generation_" + stage_counter, sub_Generation-1 , envir = .GlobalEnv) # save the Generation's number to Global Env
	# These are the information for the final expanded network at end of each stage
	# find the name of the reactions in new generated expanded network
	assign("Core_Reactions_name",Reconstructed_Reaction[Core_Reactions], envir = .GlobalEnv) # save to the Global Env
	# find the name of the metabolites in new generated expanded network
	assign("Core_Metabolites_name",Metabolites[Core_Metabolites], envir = .GlobalEnv) # save to the Global Env
	# extracting the flux for each reaction from the solutions
	assign("sol_obj",last_obj, envir = .GlobalEnv) # save to the Global Env
	assign("Rxnsize",ncol(Smat_Core), envir = .GlobalEnv) # number of reactions in new s mat
	assign("Metsize",nrow(Smat_Core), envir = .GlobalEnv) # number of metabolites in new s mat
	assign("Generation",Generation, envir = .GlobalEnv) # save the Generation's number to Global Env
	Core_Reactions<-which(Reactions %in% Core_Reactions_name ) #return back the index to the original one
	print(c("sol_obj_" + stage_counter,   get("sol_obj_"+ stage_counter , envir = .GlobalEnv)),sep="\n",quote=FALSE);
	print(c("Rxnsize_" + stage_counter,     get("Rxnsize_" + stage_counter, envir = .GlobalEnv)),sep="\n",quote=FALSE);
	print(c("Metsize_" + stage_counter,     get("Metsize_" + stage_counter, envir = .GlobalEnv)),sep="\n",quote=FALSE);
	print(c("Sub_Generation_" + stage_counter,     get("Sub_Generation_" + stage_counter, envir = .GlobalEnv)),sep="\n",quote=FALSE);
	print(c("End of Stage : ",stage_counter  ),sep="\n",quote=FALSE);
}


#' ################################
#' ###                          ###
#' ###     Printing outputs     ###
#' ###                          ###
#' ################################


if(abs(sol_obj) > 0) {
	if(is.numeric(sol_x[1])){
		sol_x[1]<-abs(sol_x[1])
	}
	Lst_1<-list()
	Lst_2<-list()
	for(i1 in 1:ncol(Smat_Core) ) {
		Lst_1[[i1]] <-Core_Metabolites_name[which(Smat_Core[,i1 ] > 0)] # finding the production
		Lst_2[[i1]] <-Core_Metabolites_name[which(Smat_Core[,i1 ] < 0)]
	}
	# finding the consumption
	Lst_1<-unique(unlist(Lst_1))
	Lst_2<-unique(unlist(Lst_2))
	# assigning a header to the list of products
	vect_1[1:(1+length(Lst_1)) ,1 ] <-c("Production" ,Lst_1)
	# assigning a header to the list of consumption
	vect_1[1:(1+length(Lst_2)) ,2 ] <-c("Consumption",Lst_2)
	
}

Steps_Info<-do.call(rbind, List_1)

Create_Output("Production_Consumption.csv",vect_1)
Create_Output("Info_EX.csv",Steps_Info)
Create_Output("Metabolites_Generation.csv",Metabolites_Generation)
Create_Output("Reactions_Flux_Generation.csv",Reactions_Flux_Generation)
Create_Output("Reactions_Flux_Flux_RealFormula.csv",Reactions_Flux_RealFormula)

# sink(file=NULL,append=F)


