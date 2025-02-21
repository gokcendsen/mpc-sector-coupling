using JuMP
using TickTock
using XLSX
using Couenne_jll, AmplNLWriter


# ========================================================================
# Title: A cost-optimal predictive operation scheme for sector-coupled 
# energy plants with start-up delays and start-up costs
# Author: G. D. Sen, J. E. Machado, J. Schiffer
# Affiliation: [Brandenburg University of Technology Cottbus-Senftenberg]
# Description: This code implements the model predictive control (MPC) 
#              optimization scheme for sector-coupled energy plants
#              with start-up delays.
# Reference: This code is associated with the paper submitted to Smart 
#            Energy Systems for efficient and sustainable smart grids 
#            and smart cities - SENSYS [2025].
# License: MIT License 
# ========================================================================



# ESU function
# dynamical equations for electrical storage unit
function ESU(W_e, e_c, e_d, p_c, p_d)
    return W_e*(1-0.01) + (e_c*p_c - p_d/e_d)
end


# TSU function
# dynamical equations for thermal storage unit
function TSU(W_t, e_ch, e_dch, q_ch, q_dch)
    return W_t*(1-0.01) + (e_ch*q_ch - q_dch/e_dch)
end

# MPC parameters
h = 1 # time step (1 hour)
N = Int(12 / h)  # prediction horizon
time = 0:h:(48 + N * h) # simulation time (48 hours)

# Load data from Excel
data=XLSX.readdata("/SENSYS_energy_data.xlsx", "January!A2:F62")
Q_dem = data[:, 2] # heat demand
P_dem = data[:, 3] # electricity demand
P_solar = data[:, 4] # solar energy 
c_el = data[:, 5]  # dynamic electricity price
c_gas = data[:, 6]  # constant gas price


# CHP parameters
eta_chp_el, eta_chp_he = 0.9, 0.9 # chp efficiencies for electricity and heat
P_chp_min, P_chp_max = 50, 150 # chp min and max electrical output limits [kW]

# Gas boiler parameters
eta_gb = 0.9 # gas boiler efficiency
Q_gb_min, Q_gb_max = 50, 250 # gas boiler min and max output limits [kW]

# hp parameters
cop_hp = 3 # coefficient of performance of heat pump 
Q_hp_min, Q_hp_max = 50, 250 # heat pump min and max output limits [kW]

# ESU parameters
eta_c, eta_d = 0.95, 0.95  # electrical storage unit charging and discharging efficiencies
P_c_min, P_c_max = 0, 200 # min and max charging limits [kW]
P_d_min, P_d_max = 0, 200 # min and max discharging limits [kW]
W_esu_min, W_esu_max = 40.0, 400.0 # min and max capacity of electrical storage unit [kW]
W_esu0 = 60.0  # initial amount of energy in ESU [kW]

# TSU parameters
eta_ch, eta_dch = 0.95, 0.95 # thermal storage unit charging and discharging efficiencies
Q_ch_min, Q_ch_max = 0, 200 # min and max charging limits [kW]
Q_dch_min, Q_dch_max = 0, 200 # min and max discharging limits [kW]
W_tsu_min, W_tsu_max = 60.0, 600.0  # min and max capacity of thermal storage unit [kW]
W_tsu0 = 80.0   # initial amount of energy in TSU [kW]

# Purchased electricity limits
p_buy_min, p_buy_max = 0, 400 # min and max electricity purchase limits [kW]


# Delay parameters
# To avoid numerical errors these parameters are not given as integer, for example Delay_chp = 1.5
# so that the CHP unit has to wait two steps and then it can produce nonzero output.
Delay_chp, Delay_gb = 1.5, 2.5 # delays [h]
UT_chp, UT_gb = 2.5, 3.5 # on time limits [h]
DT_chp, DT_gb = 0.5, 1.5 # off time limits [h]

# Hot, Warm, Cold Start time
HS_time, WS_time = 2.5, 7.5


# Closed-loop variables (arrays for storing closed-loop simulation data)
P_chp_cl, Q_chp_cl, Q_gb_cl = [], [], []
P_in_hp_cl, P_c_cl, P_d_cl, Q_ch_cl, Q_dch_cl = [], [], [], [], []
W_esu_cl, W_tsu_cl = [W_esu0], [W_tsu0]
P_buy_cl, Jcost_cl = [], []
SU_chp_cl, SU_hp_cl, SU_gb_cl = [], [], []
SD_chp_cl, SD_hp_cl, SD_gb_cl = [], [], []
V_chp_cl, V_hp_cl, V_gb_cl = [], [], []
X_ON_chp_cl, X_ON_hp_cl, X_ON_gb_cl = [], [], []
X_OFF_chp_cl, X_OFF_hp_cl, X_OFF_gb_cl = [], [], []
D_SU_chp_cl, D_SU_gb_cl = [], []
HS_chp_cl, HS_hp_cl, HS_gb_cl = [], [], []
WS_chp_cl, WS_hp_cl, WS_gb_cl = [], [], []
CS_chp_cl, CS_hp_cl, CS_gb_cl = [], [], []



# Start simulation timer 
tick()
mpciter = 0 # MPC iteration step
while mpciter <= length(time) - N - 2 
    
    model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe)) #Couenne; an MINLP solver
    
    # defining the optimization variables
    @variable(model, Q_gb[1:N] >= 0) # heat output of the gas boiler
    @variable(model, P_chp[1:N]>= 0) # electrical output of the CHP
    @variable(model, Q_chp[1:N] >= 0) # heat output of the CHP
    @variable(model, P_in_hp[1:N]>= 0) # electrical input of the heat pump
    @variable(model, P_c[1:N]>= 0) # charging variable of the ESU
    @variable(model, P_d[1:N]>= 0) # discharging variable of the ESU
    @variable(model, Q_ch[1:N]>= 0) # charging variable of the TSU
    @variable(model, Q_dch[1:N]>= 0) # discharging variable of the TSU
    @variable(model, P_buy[1:N]>= 0) # electricity purchase variable
    @variable(model, V_chp[1:N], Bin) # Boolean variable for the status of the CHP
    @variable(model, V_hp[1:N], Bin) # Boolean variable for the status of the heat pump
    @variable(model, V_gb[1:N], Bin) # Boolean variable for the status of the gas boiler
    @variable(model, D_SU_chp[1:N], Bin) # Boolean variable for start-up delay of the CHP
    @variable(model, D_SU_gb[1:N], Bin) # Boolean variable for start-up delay of the gas boiler
    @variable(model, HS_chp[1:N], Bin) # Boolean variable for hot start-up cost of the CHP
    @variable(model, HS_hp[1:N], Bin) # Boolean variable for hot start-up cost of the heat pump
    @variable(model, HS_gb[1:N], Bin) # Boolean variable for hot start-up cost of the gas boiler
    @variable(model, WS_chp[1:N], Bin) # Boolean variable for warm start-up cost of the CHP
    @variable(model, WS_hp[1:N], Bin) # Boolean variable for warm start-up cost of the heat pump
    @variable(model, WS_gb[1:N], Bin) # Boolean variable for warm start-up cost of the gas boiler
    @variable(model, CS_chp[1:N], Bin) # Boolean variable for cold start-up cost of the CHP
    @variable(model, CS_hp[1:N], Bin) # Boolean variable for cold start-up cost of the heat pump
    @variable(model, CS_gb[1:N], Bin) # Boolean variable for cold start-up cost of the gas boiler


    # global variables 
    global W_esu = [] # State of the ESU
    global W_tsu = [] # State of the TSU
    global SU_chp = [] # Boolean variable for start-up status of the CHP
    global SU_hp = [] # Boolean variable for start-up status of the heat pump
    global SU_gb = [] # Boolean variable for start-up status of the gas boiler
    global SD_chp = [] # Boolean variable for shut-down status of the CHP
    global SD_hp = [] # Boolean variable for shut-down status of the heat pump
    global SD_gb = [] # Boolean variable for shut-down status of the gas boiler
    global X_ON_chp = [] # Integer variable for on-time duration of the CHP
    global X_ON_hp = [] # Integer variable for on-time duration of the heat pump
    global X_ON_gb = [] # Integer variable for on-time duration of the gas boiler
    global X_OFF_chp = [] # Integer variable for off-time duration of the CHP
    global X_OFF_hp = [] # Integer variable for off-time duration of the heat pump
    global X_OFF_gb = [] # Integer variable for off-time duration of the gas boiler

    Jcost = 0.0 # initializing the cost value

    for i in 1:N # for loop for MPC (the problem is solved for N steps of prediction horizon) 
        
        # Dynamic constraints for storage units
        if i == 1 
            global W_esu = push!(W_esu, ESU(W_esu0, eta_c, eta_d, P_c[i], P_d[i]))
            global W_tsu = push!(W_tsu, TSU(W_tsu0, eta_ch, eta_dch, Q_ch[i], Q_dch[i]))
        else
            global W_esu = push!(W_esu, ESU(W_esu[i-1], eta_c, eta_d, P_c[i], P_d[i]))
            global W_tsu = push!(W_tsu, TSU(W_tsu[i-1], eta_ch, eta_dch, Q_ch[i], Q_dch[i]))
        end
        
        
        # Storage units cannot charge and discharge simultaneously. This constraint is to prevent that
        @constraint(model, P_c[i] * P_d[i] <= 0)
        @constraint(model, (Q_ch[i]) * Q_dch[i] <= 0)
        
        # Electrical power balance constraint
        @constraint(model, P_solar[mpciter + i] + P_chp[i] + P_d[i] + P_buy[i] - P_c[i] - P_in_hp[i] - P_dem[mpciter + i]  == 0)

        # Thermal power balance constraint
        @constraint(model,  Q_chp[i] + Q_gb[i] + cop_hp * P_in_hp[i] + Q_dch[i] - Q_ch[i]- Q_dem[mpciter + i] == 0 )
  
        # Add start-up and operating constraints based on iteration (mpciter) logic
        # Constraints for the first step of the MPC iteration
        if mpciter == 0 && i == 1

            #ramp up/down constraint for the CHP unit
            @constraint(model, 0 <= P_chp[i] <= 50)
            @constraint(model, 0 <= Q_gb[i] <= 100)

            #Start Up booleans
            global SU_chp = push!(SU_chp, V_chp[i])
            global SU_hp = push!(SU_hp, V_hp[i])
            global SU_gb = push!(SU_gb, V_gb[i])

            #Shut Down booleans
            global SD_chp = push!(SD_chp, 0)
            global SD_hp = push!(SD_hp, 0)
            global SD_gb = push!(SD_gb, 0)

            #ON/OFF time duration
            global X_ON_chp = push!(X_ON_chp, V_chp[i])
            global X_ON_hp = push!(X_ON_hp, V_hp[i])
            global X_ON_gb = push!(X_ON_gb, V_gb[i])

            global X_OFF_chp = push!(X_OFF_chp, (2+1)*(1-V_chp[i])) # we assume that all the units were off for 2 hours at the begining of the simulation.
            global X_OFF_hp = push!(X_OFF_hp, (2+1)*(1-V_hp[i]))
            global X_OFF_gb = push!(X_OFF_gb, (2+1)*(1-V_gb[i]))

            # Hot Start, Warm Start, Cold Start Booleans
            @constraint(model, WS_chp[i] == 0)
            @constraint(model, CS_chp[i] == 0)
            @constraint(model, HS_chp[i] + WS_chp[i] + CS_chp[i] == SU_chp[i])

            @constraint(model, WS_hp[i] == 0)
            @constraint(model, CS_hp[i] == 0)
            @constraint(model, HS_hp[i] + WS_hp[i] + CS_hp[i] == SU_hp[i])

            @constraint(model, WS_gb[i] == 0)
            @constraint(model, CS_gb[i] == 0)
            @constraint(model, HS_gb[i] + WS_gb[i] + CS_gb[i] == SU_gb[i])

            # Start Up delay constraints 
            # In the first step of the MPC iteration, they are zero. 
            @constraint(model, D_SU_chp[i] <= 0)
            @constraint(model, D_SU_gb[i] <= 0)
        
        elseif mpciter == 0 && i > 1
            #ramp up/down constraint for the CHP unit
            @constraint(model, -50 <= P_chp[i] - P_chp[i-1] <= 50)
            @constraint(model, -100 <= Q_gb[i] - Q_gb[i-1] <= 100)
            #Start Up booleans
            global SU_chp = push!(SU_chp, V_chp[i] * (1 - V_chp[i-1]))
            global SU_hp = push!(SU_hp, V_hp[i] * (1 - V_hp[i-1]))
            global SU_gb = push!(SU_gb, V_gb[i] * (1 - V_gb[i-1]))

            #Shut Down booleans
            global SD_chp = push!(SD_chp,  V_chp[i-1] * (1 - V_chp[i]))
            global SD_hp = push!(SD_hp, V_hp[i-1] * (1 - V_hp[i]))
            global SD_gb = push!(SD_gb, V_gb[i-1] * (1 - V_gb[i]))


            #ON/OFF time duration
            global X_ON_chp = push!(X_ON_chp,  (X_ON_chp[i-1]+1)*V_chp[i-1]*V_chp[i]+V_chp[i]*(1-V_chp[i-1]))
            global X_ON_hp = push!(X_ON_hp, (X_ON_hp[i-1]+1)*V_hp[i-1]*V_hp[i]+V_hp[i]*(1-V_hp[i-1]))
            global X_ON_gb = push!(X_ON_gb, (X_ON_gb[i-1]+1)*V_gb[i-1]*V_gb[i]+V_gb[i]*(1-V_gb[i-1]))

            global X_OFF_chp = push!(X_OFF_chp, (X_OFF_chp[i-1]+1)*(1-V_chp[i-1])*(1-V_chp[i])+V_chp[i-1]*(1-V_chp[i]))
            global X_OFF_hp = push!(X_OFF_hp, (X_OFF_hp[i-1]+1)*(1-V_hp[i-1])*(1-V_hp[i])+V_hp[i-1]*(1-V_hp[i]))
            global X_OFF_gb = push!(X_OFF_gb, (X_OFF_gb[i-1]+1)*(1-V_gb[i-1])*(1-V_gb[i])+V_gb[i-1]*(1-V_gb[i]))

            #Minimum UP/DOWN time
            @constraint(model, (X_ON_chp[i-1]-UT_chp)*(V_chp[i-1]-V_chp[i])>=0)
            @constraint(model, (X_ON_gb[i-1]-UT_gb)*(V_gb[i-1]-V_gb[i])>=0)

            @constraint(model, (X_OFF_chp[i-1]-DT_chp)*(V_chp[i]-V_chp[i-1])>=0)
            @constraint(model, (X_OFF_gb[i-1]-DT_gb)*(V_gb[i]-V_gb[i-1])>=0)

            # Hot Start, Warm Start, Cold Start Booleans
            @constraint(model, (HS_time - X_OFF_chp[i-1])*HS_chp[i] >= 0)
            @constraint(model, (X_OFF_chp[i-1] - HS_time)*(WS_time-X_OFF_chp[i-1])*WS_chp[i] >= 0)
            @constraint(model, (X_OFF_chp[i-1] - WS_time)*CS_chp[i] >= 0)
            @constraint(model, HS_chp[i] + WS_chp[i] + CS_chp[i] == SU_chp[i])

            @constraint(model, (HS_time - X_OFF_hp[i-1])*HS_hp[i] >= 0)
            @constraint(model, (X_OFF_hp[i-1] - HS_time)*(WS_time-X_OFF_hp[i-1])*WS_hp[i] >= 0)
            @constraint(model, (X_OFF_hp[i-1] - WS_time)*CS_hp[i] >= 0)
            @constraint(model, HS_hp[i] + WS_hp[i] + CS_hp[i] == SU_hp[i])

            @constraint(model, (HS_time - X_OFF_gb[i-1])*HS_gb[i] >= 0)
            @constraint(model, (X_OFF_gb[i-1] - HS_time)*(WS_time-X_OFF_gb[i-1])*WS_gb[i] >= 0)
            @constraint(model, (X_OFF_gb[i-1] - WS_time)*CS_gb[i] >= 0)
            @constraint(model, HS_gb[i] + WS_gb[i] + CS_gb[i] == SU_gb[i])
       

            # Start Up delay constraints  
            @constraint(model, (Delay_chp-X_ON_chp[i-1])*(V_chp[i]-D_SU_chp[i])>=0)
            @constraint(model, (X_ON_chp[i-1]-Delay_chp)*(D_SU_chp[i])>=0)
            @constraint(model, (Delay_gb-X_ON_gb[i-1])*(V_gb[i]-D_SU_gb[i]) >=0)
            @constraint(model, (X_ON_gb[i-1]-Delay_gb)*(D_SU_gb[i])>=0)
            @constraint(model, D_SU_chp[i]<=V_chp[i])
            @constraint(model, D_SU_gb[i]<=V_gb[i])

    
        # Constraints for the following steps
        elseif mpciter > 0 && i == 1

            #ramp up/down constraint for the CHP unit
            @constraint(model, -50 <= P_chp[i] - P_chp_cl[mpciter] <= 50)
            @constraint(model, -100 <= Q_gb[i] -Q_gb_cl[mpciter] <= 100)

            #Start Up booleans
            global SU_chp = push!(SU_chp, V_chp[i] * (1 - V_chp_cl[mpciter]))
            global SU_hp = push!(SU_hp, V_hp[i] * (1 - V_hp_cl[mpciter]))
            global SU_gb = push!(SU_gb, V_gb[i] * (1 - V_gb_cl[mpciter]))

            #Shut Down booleans
            global SD_chp = push!(SD_chp,  V_chp_cl[mpciter] * (1 - V_chp[i]))
            global SD_hp = push!(SD_hp, V_hp_cl[mpciter] * (1 - V_hp[i]))
            global SD_gb = push!(SD_gb, V_gb_cl[mpciter] * (1 - V_gb[i]))

            #ON/OFF time duration
            global X_ON_chp = push!(X_ON_chp,  (X_ON_chp_cl[mpciter]+1)*V_chp_cl[mpciter]*V_chp[i]+V_chp[i]*(1-V_chp_cl[mpciter]))
            global X_ON_hp = push!(X_ON_hp, (X_ON_hp_cl[mpciter]+1)*V_hp_cl[mpciter]*V_hp[i]+V_hp[i]*(1-V_hp_cl[mpciter]))
            global X_ON_gb = push!(X_ON_gb, (X_ON_gb_cl[mpciter]+1)*V_gb_cl[mpciter]*V_gb[i]+V_gb[i]*(1-V_gb_cl[mpciter]))

            global X_OFF_chp = push!(X_OFF_chp, (X_OFF_chp_cl[mpciter]+1)*(1-V_chp_cl[mpciter])*(1-V_chp[i])+V_chp_cl[mpciter]*(1-V_chp[i]))
            global X_OFF_hp = push!(X_OFF_hp, (X_OFF_hp_cl[mpciter]+1)*(1-V_hp_cl[mpciter])*(1-V_hp[i])+V_hp_cl[mpciter]*(1-V_hp[i]))
            global X_OFF_gb = push!(X_OFF_gb, (X_OFF_gb_cl[mpciter]+1)*(1-V_gb_cl[mpciter])*(1-V_gb[i])+V_gb_cl[mpciter]*(1-V_gb[i]))

            #Minimum UP/DOWN time
            @constraint(model, (X_ON_chp_cl[mpciter]-UT_chp)*(V_chp_cl[mpciter]-V_chp[i])>=0)
            @constraint(model, (X_ON_gb_cl[mpciter]-UT_gb)*(V_gb_cl[mpciter]-V_gb[i])>=0)

            @constraint(model, (X_OFF_chp_cl[mpciter]-DT_chp)*(V_chp[i]-V_chp_cl[mpciter])>=0)
            @constraint(model, (X_OFF_gb_cl[mpciter]-DT_gb)*(V_gb[i]-V_gb_cl[mpciter])>=0)

            # Hot Start, Warm Start, Cold Start Booleans
            @constraint(model, (HS_time - X_OFF_chp_cl[mpciter])*HS_chp[i] >= 0)
            @constraint(model, (X_OFF_chp_cl[mpciter] - HS_time)*(WS_time-X_OFF_chp_cl[mpciter])*WS_chp[i] >= 0)
            @constraint(model, (X_OFF_chp_cl[mpciter] - WS_time)*CS_chp[i] >= 0)
            @constraint(model, HS_chp[i] + WS_chp[i] + CS_chp[i] == SU_chp[i])

            @constraint(model, (HS_time - X_OFF_hp_cl[mpciter])*HS_hp[i] >= 0)
            @constraint(model, (X_OFF_hp_cl[mpciter] - HS_time)*(WS_time-X_OFF_hp_cl[mpciter])*WS_hp[i] >= 0)
            @constraint(model, (X_OFF_hp_cl[mpciter] - WS_time)*CS_hp[i] >= 0)
            @constraint(model, HS_hp[i] + WS_hp[i] + CS_hp[i] == SU_hp[i])

            @constraint(model, (HS_time - X_OFF_gb_cl[mpciter])*HS_gb[i] >= 0)
            @constraint(model, (X_OFF_gb_cl[mpciter] - HS_time)*(WS_time-X_OFF_gb_cl[mpciter])*WS_gb[i] >= 0)
            @constraint(model, (X_OFF_gb_cl[mpciter] - WS_time)*CS_gb[i] >= 0)
            @constraint(model, HS_gb[i] + WS_gb[i] + CS_gb[i] == SU_gb[i])

            #Start Up delay constraints   
            @constraint(model, (Delay_chp-X_ON_chp_cl[mpciter])*(V_chp[i]-D_SU_chp[i])>=0)
            @constraint(model, (X_ON_chp_cl[mpciter]-Delay_chp)*(D_SU_chp[i])>=0)
            @constraint(model, (Delay_gb-X_ON_gb_cl[mpciter])*(V_gb[i]-D_SU_gb[i]) >=0)
            @constraint(model, (X_ON_gb_cl[mpciter]-Delay_gb)*(D_SU_gb[i])>=0)
            @constraint(model, D_SU_chp[i]<=V_chp[i])
            @constraint(model, D_SU_gb[i]<=V_gb[i])
        
        elseif mpciter > 0 && i > 1
            #ramp up/down constraint for the CHP unit
            @constraint(model, -50 <= P_chp[i] - P_chp[i-1] <= 50)
            @constraint(model, -100 <= Q_gb[i] - Q_gb[i-1] <= 100)
            #Start Up booleans
            global SU_chp = push!(SU_chp, V_chp[i] * (1 - V_chp[i-1]))
            global SU_hp = push!(SU_hp, V_hp[i] * (1 - V_hp[i-1]))
            global SU_gb = push!(SU_gb, V_gb[i] * (1 - V_gb[i-1]))

            #Shut Down booleans
            global SD_chp = push!(SD_chp,  V_chp[i-1] * (1 - V_chp[i]))
            global SD_hp = push!(SD_hp, V_hp[i-1] * (1 - V_hp[i]))
            global SD_gb = push!(SD_gb, V_gb[i-1] * (1 - V_gb[i]))

            #ON/OFF time duration
            global X_ON_chp = push!(X_ON_chp,  (X_ON_chp[i-1]+1)*V_chp[i-1]*V_chp[i]+V_chp[i]*(1-V_chp[i-1]))
            global X_ON_hp = push!(X_ON_hp, (X_ON_hp[i-1]+1)*V_hp[i-1]*V_hp[i]+V_hp[i]*(1-V_hp[i-1]))
            global X_ON_gb = push!(X_ON_gb, (X_ON_gb[i-1]+1)*V_gb[i-1]*V_gb[i]+V_gb[i]*(1-V_gb[i-1]))

            global X_OFF_chp = push!(X_OFF_chp, (X_OFF_chp[i-1]+1)*(1-V_chp[i-1])*(1-V_chp[i])+V_chp[i-1]*(1-V_chp[i]))
            global X_OFF_hp = push!(X_OFF_hp, (X_OFF_hp[i-1]+1)*(1-V_hp[i-1])*(1-V_hp[i])+V_hp[i-1]*(1-V_hp[i]))
            global X_OFF_gb = push!(X_OFF_gb, (X_OFF_gb[i-1]+1)*(1-V_gb[i-1])*(1-V_gb[i])+V_gb[i-1]*(1-V_gb[i]))

            #Minimum UP/DOWN time
            @constraint(model, (X_ON_chp[i-1]-UT_chp)*(V_chp[i-1]-V_chp[i])>=0)
            @constraint(model, (X_ON_gb[i-1]-UT_gb)*(V_gb[i-1]-V_gb[i])>=0)

            @constraint(model, (X_OFF_chp[i-1]-DT_chp)*(V_chp[i]-V_chp[i-1])>=0)
            @constraint(model, (X_OFF_gb[i-1]-DT_gb)*(V_gb[i]-V_gb[i-1])>=0)

            # Hot Start, Warm Start, Cold Start Booleans
            @constraint(model, (HS_time - X_OFF_chp[i-1])*HS_chp[i] >= 0)
            @constraint(model, (X_OFF_chp[i-1] - HS_time)*(WS_time-X_OFF_chp[i-1])*WS_chp[i] >= 0)
            @constraint(model, (X_OFF_chp[i-1] - WS_time)*CS_chp[i] >= 0)
            @constraint(model, HS_chp[i] + WS_chp[i] + CS_chp[i] == SU_chp[i])

            @constraint(model, (HS_time - X_OFF_hp[i-1])*HS_hp[i] >= 0)
            @constraint(model, (X_OFF_hp[i-1] - HS_time)*(WS_time-X_OFF_hp[i-1])*WS_hp[i] >= 0)
            @constraint(model, (X_OFF_hp[i-1] - WS_time)*CS_hp[i] >= 0)
            @constraint(model, HS_hp[i] + WS_hp[i] + CS_hp[i] == SU_hp[i])

            @constraint(model, (HS_time - X_OFF_gb[i-1])*HS_gb[i] >= 0)
            @constraint(model, (X_OFF_gb[i-1] - HS_time)*(WS_time-X_OFF_gb[i-1])*WS_gb[i] >= 0)
            @constraint(model, (X_OFF_gb[i-1] - WS_time)*CS_gb[i] >= 0)
            @constraint(model, HS_gb[i] + WS_gb[i] + CS_gb[i] == SU_gb[i])

            #Start Up delay constraints    
            @constraint(model, (Delay_chp-X_ON_chp[i-1])*(V_chp[i]-D_SU_chp[i])>=0)
            @constraint(model, (X_ON_chp[i-1]-Delay_chp)*(D_SU_chp[i])>=0)
            @constraint(model, (Delay_gb-X_ON_gb[i-1])*(V_gb[i]-D_SU_gb[i]) >=0)
            @constraint(model, (X_ON_gb[i-1]-Delay_gb)*(D_SU_gb[i])>=0)
            @constraint(model, D_SU_chp[i]<=V_chp[i])
            @constraint(model, D_SU_gb[i]<=V_gb[i])
        
        end

        # Output constraints for CHP, Heat Pump and Gas Boiler
        @constraint(model, D_SU_chp[i]*P_chp_min <= P_chp[i])
        @constraint(model, P_chp[i] <= P_chp_max*D_SU_chp[i])
        @constraint(model, D_SU_gb[i]*Q_gb_min <= Q_gb[i])
        @constraint(model, Q_gb[i] <= Q_gb_max*D_SU_gb[i])
        @constraint(model, V_hp[i]*Q_hp_min <= cop_hp*P_in_hp[i])
        @constraint(model, cop_hp*P_in_hp[i] <= V_hp[i]*Q_hp_max)


        # Cost function summing
        Jcost += h * c_gas[mpciter + i] * (P_chp[i] / eta_chp_el + Q_chp[i] / eta_chp_he) # gas costs that chp used
        Jcost += h * c_gas[mpciter + i] * (Q_gb[i] / eta_gb) # gas costs that gas boiler used 
        Jcost += h * c_el[mpciter + i] * P_buy[i] # purchased electricity costs
        Jcost += 15 * HS_chp[i] + 30 * WS_chp[i] +  45 * CS_chp[i]  # start up costs for chp
        Jcost += 5 * HS_hp[i]  + 10 * WS_hp[i] + 20 * CS_hp[i]  # start up costs for heat pump
        Jcost += 10 * HS_gb[i]  + 20 * WS_gb[i] + 30 * CS_gb[i]    # start up costs for gas boiler
        Jcost += 15 * SD_chp[i] # shut down cost for chp
        Jcost += 5 * SD_hp[i] # shut down cost for heat pump
        Jcost += 10 * SD_gb[i]  # shut down cost for gas boiler
    end
    
    
    # Additional constraints on variables
    @constraint(model,  1.7*P_chp == Q_chp) # heat output of the CHP unit should be 1.7 times electrical output of the CHP unit
 
    @constraint(model, W_esu_min  .<= W_esu .<=  W_esu_max) # min and max capacity constraint of ESU
    @constraint(model, W_tsu_min  .<= W_tsu .<=  W_tsu_max) # min and max capacity constraint of TSU
    @constraint(model, P_c_min .<= P_c .<= P_c_max) # min and max charging capacity constraint of ESU
    @constraint(model, P_d_min .<= P_d .<= P_d_max) # min and max discharging capacity constraint of ESU
    @constraint(model, Q_ch_min .<= Q_ch .<= Q_ch_max) # min and max charging capacity constraint of TSU
    @constraint(model, Q_dch_min .<= Q_dch .<= Q_dch_max) # min and max discharging capacity constraint of TSU 
    @constraint(model, p_buy_min .<= P_buy .<= p_buy_max) # min and max electricity purchase constraint
    
    @constraint(model, W_esu[N] == W_esu_max ) # terminal constraints for feasibility
    @constraint(model, W_tsu[N] == W_tsu_max )

    # Objective function
    @objective(model, Min, Jcost)
    optimize!(model)

 
    # Update MPC closed-loop solutions
    # all the results are stored in arrays
    # in each iteration of MPC we only apply the first element of the optimal results 
    # therefore we only store the first element 
    push!(P_chp_cl, value(P_chp[1]))
    push!(Q_chp_cl, value(Q_chp[1]))
    push!(Q_gb_cl, value(Q_gb[1]))
    push!(P_in_hp_cl, value(P_in_hp[1]))
    push!(P_c_cl, value(P_c[1]))
    push!(P_d_cl, value(P_d[1]))
    push!(Q_ch_cl, value(Q_ch[1]))  
    push!(Q_dch_cl, value(Q_dch[1]))    
    push!(P_buy_cl, value(P_buy[1]))


    if mpciter == 0
        global SU_chp_cl = push!(SU_chp_cl, round.(value(V_chp[1])))
        global SU_hp_cl = push!(SU_hp_cl, round.(value(V_hp[1])))
        global SU_gb_cl = push!(SU_gb_cl, round.(value(V_gb[1])))

        global SD_chp_cl = push!(SD_chp_cl,0)
        global SD_hp_cl = push!(SD_hp_cl, 0)
        global SD_gb_cl = push!(SD_gb_cl, 0)

        #ON/OFF time duration
        global X_ON_chp_cl= push!(X_ON_chp_cl, round.(value(V_chp[1])))
        global X_ON_hp_cl = push!(X_ON_hp_cl, round.(value(V_hp[1])))
        global X_ON_gb_cl = push!(X_ON_gb_cl, round.(value(V_gb[1])))

        global X_OFF_chp_cl = push!(X_OFF_chp_cl, round.((2+1)*(1-value(V_chp[1]))))
        global X_OFF_hp_cl = push!(X_OFF_hp_cl, round.((2+1)*(1-value(V_hp[1]))))
        global X_OFF_gb_cl = push!(X_OFF_gb_cl, round.((2+1)*(1-value(V_gb[1]))))
    else
        global SU_chp_cl = push!(SU_chp_cl, round.(value(V_chp[1]) * (1 - V_chp_cl[mpciter])))
        global SU_hp_cl = push!(SU_hp_cl, round.(value(V_hp[1]) * (1 - V_hp_cl[mpciter])))
        global SU_gb_cl = push!(SU_gb_cl, round.(value(V_gb[1]) * (1 - V_gb_cl[mpciter])))

        global SD_chp_cl = push!(SD_chp_cl, round.(V_chp_cl[mpciter]* (1 - value(V_chp[1]))))
        global SD_hp_cl = push!(SD_hp_cl, round.(V_hp_cl[mpciter]* (1 - value(V_hp[1]))))
        global SD_gb_cl = push!(SD_gb_cl, round.(V_gb_cl[mpciter]* (1 - value(V_gb[1]))))

        #ON/OFF time duration
        global X_ON_chp_cl = push!(X_ON_chp_cl, round.((X_ON_chp_cl[mpciter]+1)*V_chp_cl[mpciter]*value(V_chp[1])+value(V_chp[1])*(1-V_chp_cl[mpciter])))
        global X_ON_hp_cl = push!(X_ON_hp_cl,round.((X_ON_hp_cl[mpciter]+1)*V_hp_cl[mpciter]*value(V_hp[1])+value(V_hp[1])*(1-V_hp_cl[mpciter])))
        global X_ON_gb_cl = push!(X_ON_gb_cl,round.((X_ON_gb_cl[mpciter]+1)*V_gb_cl[mpciter]*value(V_gb[1])+value(V_gb[1])*(1-V_gb_cl[mpciter])))

        global X_OFF_chp_cl = push!(X_OFF_chp_cl, round.((X_OFF_chp_cl[mpciter]+1)*(1-V_chp_cl[mpciter])*(1-value(V_chp[1]))+V_chp_cl[mpciter]*(1-value(V_chp[1]))))
        global X_OFF_hp_cl = push!(X_OFF_hp_cl, round.((X_OFF_hp_cl[mpciter]+1)*(1-V_hp_cl[mpciter])*(1-value(V_hp[1]))+V_hp_cl[mpciter]*(1-value(V_hp[1]))))
        global X_OFF_gb_cl = push!(X_OFF_gb_cl, round.((X_OFF_gb_cl[mpciter]+1)*(1-V_gb_cl[mpciter])*(1-value(V_gb[1]))+V_gb_cl[mpciter]*(1-value(V_gb[1]))))
    end


    global V_chp_cl = push!(V_chp_cl, round.(value(V_chp[1])))
    global V_hp_cl = push!(V_hp_cl, round.(value(V_hp[1])))
    global V_gb_cl = push!(V_gb_cl, round.(value(V_gb[1])))
    

    push!(HS_chp_cl, round.(value(HS_chp[1])))
    push!(HS_hp_cl, round.(value(HS_hp[1])))
    push!(HS_gb_cl, round.(value(HS_gb[1])))

    push!(WS_chp_cl, round.(value(WS_chp[1])))
    push!(WS_hp_cl, round.(value(WS_hp[1])))
    push!(WS_gb_cl, round.(value(WS_gb[1])))

    push!(CS_chp_cl, round.(value(CS_chp[1])))
    push!(CS_hp_cl, round.(value(CS_hp[1])))
    push!(CS_gb_cl, round.(value(CS_gb[1])))

    push!(Jcost_cl, objective_value(model))
    push!(D_SU_chp_cl, round.(value(D_SU_chp[1])))
    push!(D_SU_gb_cl, round.(value(D_SU_gb[1])))


     # Update stored energy values
     global W_esu0 = ESU(W_esu0, eta_c, eta_d, value(P_c[1]), value(P_d[1]))
     global W_tsu0 = TSU(W_tsu0, eta_ch, eta_dch, value(Q_ch[1]), value(Q_dch[1]))

     push!(W_esu_cl, W_esu0)
     push!(W_tsu_cl, W_tsu0)
 
     # Increment MPC iteration
     global mpciter += 1
     show(mpciter)
end
tock()
