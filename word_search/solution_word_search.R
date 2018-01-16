# Library calls
library(data.table)

## DATA FOR TESTING IN testing_data.Rdata

## STUDENT SOL CODE

find_firstchar <- function(mtx, first_char) {
  return(which(mtx == first_char, arr.ind = TRUE))
}

make_poss_coords <- function(mtx, prev_char_pos) { # 3.7
  possCoords <- data.frame(rows = c(prev_char_pos[1], prev_char_pos[1] - 1, prev_char_pos[1] + 1),
                           cols = c(prev_char_pos[2], prev_char_pos[2] - 1, prev_char_pos[2] + 1))

  combns <- lapply(possCoords$rows, function(x) data.frame(cbind(x, possCoords$cols)))
  combns <- rbindlist(combns)

  #Take out the base coordinate - don't want to end up going backwarmx
  rmv_prev_char_pos <- apply(combns, 1, function(x) x != prev_char_pos)

  filter_prev_char_pos <- apply(rmv_prev_char_pos, 2, function(x) !all(x == FALSE))

  searchCoords <- combns[filter_prev_char_pos, ]

  #Take out coordinates not in the data frame
  names(searchCoords) <- c("rows", "cols")
  rowLimsCheck <- searchCoords$rows %in% seq(1, nrow(mtx))
  colLimsCheck <- searchCoords$cols %in% seq(1, ncol(mtx))

  searchCoords <- searchCoords[rowLimsCheck & colLimsCheck, ]

  return(searchCoords)
}

eval_coords <- function(mtx, all_poss_coords, prev1_char_pos, prev2_char_pos, next_char) { #3.7
  if (!all(is.na(prev2_char_pos))) {  # if both previous chars are defined
    direct <- prev1_char_pos - prev2_char_pos
    poss_next_char_pos <- prev1_char_pos + direct
    next_char_pos <- all_poss_coords[unlist(all_poss_coords[,1])==unlist(poss_next_char_pos[1]) &
                                         unlist(all_poss_coords[,2])==unlist(poss_next_char_pos[2]), ]
    if(nrow(next_char_pos)>0 && next_char == mtx[next_char_pos[[1]],next_char_pos[[2]]]) {
      return(next_char_pos)
    } else {
      return(NA)
    }
  } else { # if only prev1 defined, find all positions for next letter
    all_poss_coords_wl <- cbind(all_poss_coords, diag(mtx[unlist(all_poss_coords[,1]),
                                                          unlist(all_poss_coords[,2])]))
    next_char_pos <- all_poss_coords_wl[unlist(all_poss_coords_wl[,3])==next_char,1:2]
    if(nrow(next_char_pos)>0) {
     # class(next_char_pos) <- "numeric"
      return(next_char_pos)
    } else {
      return(NA)
    }
  }
}

det_which_dir <- function(first_coord, second_coord) { #3.8
  ref_lib_wrd_y <- data.frame(num=c(1,0,-1),wrd=c("down","","up"))
  ref_lib_wrd_x <- data.frame(num=c(1,0,-1),wrd=c("right","","left"))
  diff <- second_coord - first_coord
  if (any(abs(diff)>1)) {
    stop("Coordinates too far apart!")
  }
  if (all(diff!=0)) {
    return(paste0(ref_lib_wrd_y[diff[1]==ref_lib_wrd_y[,1],2], "-",
                  ref_lib_wrd_x[diff[2]==ref_lib_wrd_x[,1],2]))
  } else {
     both_vec <- c(paste0(ref_lib_wrd_y[diff[1]==ref_lib_wrd_y[,1],2]),
                   paste0(ref_lib_wrd_x[diff[2]==ref_lib_wrd_x[,1],2]))
     return(both_vec[diff!=0])
  }
}

word_Search <- function(mtx, wrd) {  # 3.9
  split <- sapply(seq(1, nchar(wrd)), function(i) substring(wrd, i, i))
  first_char_locs <- find_firstchar(mtx, split[1])

  if (nrow(first_char_locs)==0) {
    stop("Word not in matrix!")
  } else {
    poss_coords_lst <- apply(first_char_locs, 1, FUN=make_poss_coords, mtx=mtx)

    res_bin <-  lapply(poss_coords_lst, function(x) list(poss_coords=x, prev1_loc=NA,
                                                       prev2_loc=NA, next_char=NA,
                                                       found_coords=NA, dir=NA))

    for (jl in 1:nrow(first_char_locs)) {
      res_bin[[jl]]$prev1_loc <- first_char_locs[jl,]
      res_bin[[jl]]$found_coords <- first_char_locs[jl,]
    }

    for (js in 2:length(split)) {

      for (jl in 1:length(res_bin)) {
        res_bin[[jl]]$next_char <- split[js]
      }

      matched_coords_lst <- lapply(res_bin, function(l) eval_coords(mtx=mtx,
                                                                    all_poss_coords =l$poss_coords,
                                                                    prev1_char_pos = l$prev1_loc,
                                                                    prev2_char_pos = l$prev2_loc,
                                                                    next_char = l$next_char)
                                  )

      if (sum(is.na(matched_coords_lst)) > 0) {
        res_bin <- res_bin[-which(is.na(matched_coords_lst), arr.ind = T)]
        matched_coords_lst <- matched_coords_lst[-which(is.na(matched_coords_lst), arr.ind = T)]
      }

      if (is.list(res_bin) & length(res_bin) == 0) {
        stop("Word not found in matrix!")
      }

      for (jl in 1:length(res_bin)) {
        if (nrow(matched_coords_lst[[jl]]) > 1) {
          res_bin[[jl]]$prev2_loc <- res_bin[[jl]]$prev1_loc
          first_extra <- length(res_bin) + 1
          last_extra <- length(res_bin) + nrow(matched_coords_lst[[jl]]) - 1
          res_bin[first_extra:last_extra] <- lapply(first_extra:last_extra,
                                                    function(x, r) r, r=res_bin[[jl]])
          jll_loop_ind <- c(first_extra: last_extra,jl)
          for (jll in jll_loop_ind) {
            df_index <- which(jll==jll_loop_ind, arr.ind = T)
            res_bin[[jll]]$prev1_loc <- unlist(matched_coords_lst[[jl]][df_index,])
            res_bin[[jll]]$found_coords <- rbind(res_bin[[jl]]$found_coords,
                                               unlist(matched_coords_lst[[jl]][df_index,]))
          }
        } else {
          res_bin[[jl]]$prev2_loc <- res_bin[[jl]]$prev1_loc
          res_bin[[jl]]$prev1_loc <- unlist(matched_coords_lst[[jl]])
          res_bin[[jl]]$found_coords <- rbind(res_bin[[jl]]$found_coords,
                                              unlist(matched_coords_lst[[jl]]))
        }
      }

      curr_char_locs <- matrix(NA, length(res_bin), 2)
      for (jl in 1:length(res_bin)) {
        curr_char_locs[jl,] <- unlist(res_bin[[jl]]$prev1_loc)
      }

      poss_coords_lst <- apply(curr_char_locs, 1, FUN=make_poss_coords, mtx=mtx)
      for (jl in 1:length(res_bin)) {
        res_bin[[jl]]$poss_coords <- poss_coords_lst[[jl]]
      }
    }
  }
  for (jl in 1:length(res_bin)) {
    res_bin[[jl]]$dir <- det_which_dir(res_bin[[jl]]$found_coords[1,],
                                       res_bin[[jl]]$found_coords[2,])
  }
  report <- list()
  for (jl in 1:length(res_bin)) {
    report[[jl]] <- paste0(wrd," starts at ","(",
                           as.vector(res_bin[[jl]]$found_coords[1,][1]),",",
                           as.vector(res_bin[[jl]]$found_coords[1,][2]),")",
                           " and goes ", res_bin[[jl]]$dir)
  }
  return(unlist(report))
}

lword_search <- function(wrd_list, mtx) { # 3.10
  return(lapply(wrd_list, word_Search, mtx=mtx))
}
