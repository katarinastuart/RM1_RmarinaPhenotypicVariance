setwd("C:/Users/User/Desktop/2020_Corona_Fun/Rm1_ToadRDA/Code")

####### Common Garden raised progeny

ToadTraits <- read.csv("CaneToad_Pheno_Raw.csv",stringsAsFactors=TRUE,sep=",")
str(ToadTraits)

names(ToadTraits)
[1] "SETUP.DETAILS.Treatment.Group..1...LeLd.2...LeHd.3...HeLd.4...HeHd" "Exercise"                                                          
[3] "Diet"                                                               "Population.Origin"                                                 
[5] "Clutch"                                                             "Individual.number"                                                 
[7] "SUL..mm."                                                           "start.SUL"                                                         
[9] "Growth"                                                             "Leg.length..mm."                                                   
[11] "Arm.length..mm."                                                    "Head.Width..mm."                                                   
[13] "Heart.Mass..g."                                                     "Liver.Mass..g."                                                    
[15] "Escape"                                                             "Struggle..number.of.kicks."                                        
[17] "Righting.time..s."                                                  "Righting.score"                                                    
[19] "Distance.traveled.in.1.min..cm."                                    "Reluctance.Score"                                                  
[21] "Boldness.Score"                                                     "Climbing.Trial.Time"                                               
[23] "Feeding.Trial"                                                      "Spontaneous.Acitivty...Mean.distance.per.frame..cm."               
[25] "Spontaneous.Acitivty...max.distance..over.all.frames..cm."          "Spontaneous.Acitivty...95.quantile..cm."                           
[27] "exhaustion.score"                                                   "Stamina.score"                                                     
[29] "Recovery.score"                               

Population <- ToadTraits$Population.Origin

ToadTraits$Population.Origin

Clutch <- ToadTraits$Clutch

Individual <- ToadTraits$Individual.number

Sul <- ToadTraits$SUL..mm.

Growth <- ToadTraits$Growth


leg.lm = lm(Leg.length..mm. ~ SUL..mm., data=ToadTraits)
summary(leg.lm) 
Leg_Length_SC <- resid(leg.lm)

arm.lm = lm(Arm.length..mm. ~ SUL..mm., data=ToadTraits)
summary(arm.lm) 
Arm_Length_SC <- resid(arm.lm)

head.lm = lm(Head.Width..mm. ~ SUL..mm., data=ToadTraits)
summary(head.lm) 
Head_Size_SC <- resid(head.lm)

heart.lm = lm(Heart.Mass..g. ~ SUL..mm., data=ToadTraits)
summary(heart.lm) 
Heart_Mass_SC <- resid(heart.lm)

liver.lm = lm(Liver.Mass..g. ~ SUL..mm., data=ToadTraits)
summary(liver.lm) 
Liver_Mass_SC <- resid(liver.lm)

escpae.lm = lm(Escape ~ SUL..mm., data=ToadTraits)
summary(escpae.lm) 
Escape_score_SC <- resid(escpae.lm)

struggle.lm = lm(Struggle..number.of.kicks. ~ SUL..mm., data=ToadTraits)
summary(struggle.lm) 
Struggle_Score <- ToadTraits$Struggle..number.of.kicks.

righttime.lm = lm(Righting.time..s. ~ SUL..mm., data=ToadTraits)
summary(righttime.lm) 
Righting_Time_SC <- resid(righttime.lm)

righteff.lm = lm(Righting.score ~ SUL..mm., data=ToadTraits)
summary(righteff.lm) 
Righting_Effort_SC <- resid(righteff.lm)

dist.lm = lm(Distance.traveled.in.1.min..cm. ~ SUL..mm., data=ToadTraits)
summary(dist.lm) 
Speed_SC <- resid(dist.lm)

Reluct.lm = lm(Reluctance.Score ~ SUL..mm., data=ToadTraits)
summary(Reluct.lm) 
Reluctance_Score_SC <- resid(Reluct.lm)

bold.lm = lm(Boldness.Score ~ SUL..mm., data=ToadTraits)
summary(bold.lm) 
Shyness_Score_SC <- resid(bold.lm)

climb.lm = lm(Climbing.Trial.Time ~ SUL..mm., data=ToadTraits)
summary(climb.lm) 
Climbing_Score <- ToadTraits$Climbing.Trial.Time

feed.lm = lm(Feeding.Trial ~ SUL..mm., data=ToadTraits)
summary(feed.lm) 
Feeding_Score <- ToadTraits$Feeding.Trial

SAmean.lm = lm(Spontaneous.Acitivty...Mean.distance.per.frame..cm. ~ SUL..mm., data=ToadTraits)
summary(SAmean.lm) 
SpontaneousActivity_Mean <- ToadTraits$Spontaneous.Acitivty...Mean.distance.per.frame..cm.

SAmax.lm = lm(Spontaneous.Acitivty...max.distance..over.all.frames..cm. ~ SUL..mm., data=ToadTraits)
summary(SAmax.lm) 
SpontaneousActivity_Max <- ToadTraits$Spontaneous.Acitivty...max.distance..over.all.frames..cm.

SAquart.lm = lm(Spontaneous.Acitivty...95.quantile..cm. ~ SUL..mm., data=ToadTraits)
summary(SAquart.lm) 
SpontaneousActivity_95Quantile_SC <- resid(SAquart.lm)

exhuast.lm = lm(exhaustion.score ~ SUL..mm., data=ToadTraits)
summary(exhuast.lm) 
Exhaustion_Score_SC <- resid(exhuast.lm)

stamina.lm = lm(Stamina.score ~ SUL..mm., data=ToadTraits)
summary(stamina.lm) 
Stamina_Score_SC <- resid(stamina.lm)

recovery.lm = lm(Recovery.score ~ SUL..mm., data=ToadTraits)
summary(recovery.lm) 
Recover_Score <- ToadTraits$Recovery.score


CaneToad_Pheno_SC <- cbind(Population,Clutch,Individual,Sul,Growth,Leg_Length_SC,Arm_Length_SC,Head_Size_SC,Heart_Mass_SC,
                           Liver_Mass_SC,Escape_score_SC,Struggle_Score,Righting_Time_SC,Righting_Effort_SC,Speed_SC,
                           Reluctance_Score_SC,Shyness_Score_SC,Climbing_Score,Feeding_Score,SpontaneousActivity_Mean,
                           SpontaneousActivity_Max,SpontaneousActivity_95Quantile_SC,Exhaustion_Score_SC,
                           Stamina_Score_SC,Recover_Score)
str(CaneToad_Pheno_SC)
CaneToad_Pheno_SC


write.csv(CaneToad_Pheno_SC, "CaneToad_Pheno_SC.csv", row.names = F, quote = F)

#Similarly done for wild-caught toads
