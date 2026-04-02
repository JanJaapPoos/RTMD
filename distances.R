distances <- function (shapes) {
   cat("Distances between snout and caudal fin landmarks in ",shapes$unit[1]," \n") 
   print(apply(shapes$landmarks[c("snout", "caudal_fin"), , ], 3, distancePointToPoint ))
}
