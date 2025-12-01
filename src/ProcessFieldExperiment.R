################################################################
# Map data from field experiment                               #
#                                                              #
################################################################

# load data
# input folder
input_folder <- "input/FieldExperiment/"
# List all Excel files in the input folder
excel_files <- list.files(input_folder, pattern = "\\.xlsx$", full.names = TRUE)
data_field_study <- NULL;
for(input_file in excel_files) {
  
  SheetNames<- excel_sheets(input_file)
 #farm and groep name
  farm_group_name <- str_split(input_file, "/")[[1]][length(str_split(input_file, "/")[[1]])] %>%
    str_remove(".xlsx")
  for(sheet_name in SheetNames) {
    if(sheet_name != "Codes"){
      message(sheet_name)
    #extract time point from sheet name
    time_point <- as.numeric(str_split(sheet_name,"SD")[[1]][2])
    data_sheet_read <- read_excel(input_file, 
                                sheet = sheet_name)
    #check if the sheet has the expected columns
    if(!all(c("Inzendnummer", "Parameter", "Methode", "Uitslag") %in% colnames(data_sheet_read))){
      #colnames are first row 
      data_sheet_read <- rbind(colnames(data_sheet_read), 
                               data_sheet_read)
      #set colnames
      colnames(data_sheet_read) <-c("Inzendnummer", "Parameter", "Methode", "Uitslag")
    }
    #add time point to data
   data_sheet_read <- data_sheet_read %>%
   select("Inzendnummer","Parameter","Methode","Uitslag")%>%
    mutate(time_point = time_point)  # remove empty columns if they exist
  #add to output
   data_field_study <- rbind(data_field_study, 
                                    cbind(data_sheet_read, farm_group_name = farm_group_name))
    }
  }
}

# remove non-numeric Uitslag values
data_field_study_numeric <- data_field_study%>%
  separate(farm_group_name, into = c("id", "Bedrijf_str", "Farm","Groep_str" ,"Group","Treatment"), sep = " ", extra = "merge")%>%
  select(-id, -Bedrijf_str, -Groep_str) %>%
  mutate(Uitslag = ifelse(Uitslag=="<1",0,Uitslag))%>%
  filter(!is.na(Uitslag)&Uitslag!="ASP"& Uitslag!="BIN"& Uitslag!="Aspecifiek"&Uitslag!="Bloed is nagestold" & Uitslag !="Te weinig serum" & Uitslag!="TWS") %>%
  mutate(Uitslag = as.numeric(Uitslag))

#aggregate to percentage below cutt-off 
data_field_study_cut_offs <- data_field_study_numeric %>%
  group_by(time_point, Parameter,Farm,Group,Treatment) %>%
  summarise(O3 = sum(Uitslag >=3, na.rm = TRUE) / n(),
            O5 = sum(Uitslag >=5, na.rm = TRUE) / n(),
            O6 = sum(Uitslag >= 6, na.rm = TRUE) / n(),
            O7 =  sum(Uitslag >= 7, na.rm = TRUE) / n(),
            .groups = 'drop')%>%
  mutate(Treatment = str_replace_all(Treatment, pattern="BI\\+booster", replacement = "BIBOOSTER"))

#re-code the treatment names
data_field_study_cut_offs <- data_field_study_cut_offs%>%
  mutate(Treatment = str_replace(Treatment, pattern = "BI\\+booster", replacement = "BIBOOSTER") )
data_field_study_cut_offs <- data_field_study_cut_offs%>%
  mutate(Treatment = str_replace(Treatment, pattern = "[ ]", replacement = "_") )
data_field_study_cut_offs <- data_field_study_cut_offs%>%
  mutate(Treatment = str_replace(Treatment, pattern = "Neg.", replacement = "Neg") )



#save to RDS-file
dir.create(paste0(input_folder,"/Processed_Data"), showWarnings = FALSE)
saveRDS(data_field_study_cut_offs,paste0(input_folder,"/Processed_Data/data_field_study_cut_offs.RDS"))

