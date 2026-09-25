
# ============================================================
# FUNÇÃO 1
# Decodificar RLE COCO não comprimido
# ============================================================

decode_coco_rle <- function(seg) {
  
  if (is.null(seg$counts) || is.null(seg$size)) {
    stop("Segmentação sem counts ou size.")
  }
  
  if (is.character(seg$counts)) {
    stop(
      "RLE comprimido encontrado. ",
      "Este código foi feito para o RLE numérico usado pelos seus JSONs."
    )
  }
  
  counts <- as.numeric(unlist(seg$counts))
  size   <- as.integer(unlist(seg$size))
  
  if (length(size) != 2L) {
    stop("segmentation$size deveria ter comprimento 2.")
  }
  
  altura  <- size[1]
  largura <- size[2]
  
  if (
    any(!is.finite(counts)) ||
    any(counts < 0) ||
    any(counts != floor(counts))
  ) {
    stop("RLE inválido.")
  }
  
  if (sum(counts) != altura * largura) {
    stop(
      "RLE incompatível com as dimensões: ",
      "sum(counts) = ", sum(counts),
      "; esperado = ", altura * largura
    )
  }
  
  # COCO RLE:
  # começa em fundo (0) e alterna 0/1
  valores <- (seq_along(counts) - 1L) %% 2L
  
  # COCO usa ordem column-major
  mascara <- matrix(
    rep(valores, times = counts),
    nrow = altura,
    ncol = largura,
    byrow = FALSE
  )
  
  mascara
}


# ============================================================
# FUNÇÃO 2
# Extrair image_id espacial a partir do nome do JPG
# ============================================================

get_spatial_image_id <- function(image_info,
                                 json_data = NULL) {
  
  candidatos <- character(0)
  
  # images$file_name
  if (!is.null(image_info$file_name)) {
    candidatos <- c(
      candidatos,
      basename(image_info$file_name)
    )
  }
  
  # info$source$source_relative
  if (
    !is.null(json_data) &&
    !is.null(json_data$info$source$source_relative)
  ) {
    candidatos <- c(
      candidatos,
      basename(json_data$info$source$source_relative)
    )
  }
  
  candidatos <- unique(candidatos)
  
  # remover extensão
  candidatos <- tools::file_path_sans_ext(candidatos)
  
  # Exemplo:
  #
  # URB_009922_fb8a6834a6.jpg
  #
  # vira:
  #
  # URB_009922_fb8a6834a6
  
  candidatos
}


# ============================================================
# FUNÇÃO 3
# Transformar UMA máscara binária em polígono espacial
# ============================================================

mask_to_spatial_polygon <- function(
    mascara,
    grid_geometry,
    crs,
    crop_to_mask = TRUE
) {
  
  altura  <- nrow(mascara)
  largura <- ncol(mascara)
  
  pixels <- which(
    mascara == 1,
    arr.ind = TRUE
  )
  
  if (nrow(pixels) == 0) {
    return(NULL)
  }
  
  # ----------------------------------------------------------
  # Bounding box espacial da imagem
  # ----------------------------------------------------------
  
  bb <- st_bbox(grid_geometry)
  
  xmin <- as.numeric(bb["xmin"])
  xmax <- as.numeric(bb["xmax"])
  ymin <- as.numeric(bb["ymin"])
  ymax <- as.numeric(bb["ymax"])
  
  dx <- (xmax - xmin) / largura
  dy <- (ymax - ymin) / altura
  
  
  # ----------------------------------------------------------
  # Otimização:
  # usar apenas a região ocupada pela máscara
  # ----------------------------------------------------------
  
  if (crop_to_mask) {
    
    row_min <- min(pixels[, "row"])
    row_max <- max(pixels[, "row"])
    
    col_min <- min(pixels[, "col"])
    col_max <- max(pixels[, "col"])
    
    mascara_local <- mascara[
      row_min:row_max,
      col_min:col_max,
      drop = FALSE
    ]
    
    # --------------------------------------------------------
    # X
    #
    # coluna 1 começa em xmin
    # --------------------------------------------------------
    
    xmin_local <- xmin +
      (col_min - 1) * dx
    
    xmax_local <- xmin +
      col_max * dx
    
    
    # --------------------------------------------------------
    # Y
    #
    # IMPORTANTE:
    #
    # imagem:
    # origem = canto superior esquerdo
    # row aumenta para BAIXO
    #
    # espaço cartográfico:
    # Y aumenta para CIMA
    # --------------------------------------------------------
    
    ymax_local <- ymax -
      (row_min - 1) * dy
    
    ymin_local <- ymax -
      row_max * dy
    
  } else {
    
    mascara_local <- mascara
    
    xmin_local <- xmin
    xmax_local <- xmax
    ymin_local <- ymin
    ymax_local <- ymax
  }
  
  
  # ----------------------------------------------------------
  # Criar raster espacial
  # ----------------------------------------------------------
  
  r <- terra::rast(
    nrows = nrow(mascara_local),
    ncols = ncol(mascara_local),
    
    xmin = xmin_local,
    xmax = xmax_local,
    
    ymin = ymin_local,
    ymax = ymax_local,
    
    crs = st_crs(crs)$wkt
  )
  
  
  # terra percorre as células por linha,
  # da esquerda para direita e de cima para baixo.
  #
  # A matriz R precisa, portanto, ser transformada
  # para ordem row-major.
  
  valores <- as.vector(
    t(mascara_local)
  )
  
  valores[valores == 0] <- NA
  
  terra::values(r) <- valores
  
  
  # ----------------------------------------------------------
  # Raster -> vetor
  # ----------------------------------------------------------
  
  p <- terra::as.polygons(
    r,
    aggregate = TRUE,
    values = TRUE,
    na.rm = TRUE
  )
  
  if (nrow(p) == 0) {
    return(NULL)
  }
  
  p <- st_as_sf(p)
  
  # Uma anotação pode ter várias partes desconectadas.
  # Mantemos tudo como uma única MULTIPOLYGON quando necessário.
  
  geom <- st_union(
    st_geometry(p)
  )
  
  geom <- st_make_valid(geom)
  
  geom
}


# ============================================================
# FUNÇÃO 4
# Converter UM annotations.coco.json
# ============================================================

convert_one_coco <- function(
    json_file,
    grids,
    grid_index,
    clip_to_grid = FALSE
) {
  
  json_data <- jsonlite::fromJSON(
    json_file,
    simplifyVector = FALSE
  )
  
  if (
    is.null(json_data$images) ||
    length(json_data$images) == 0
  ) {
    warning("Nenhuma imagem em: ", json_file)
    return(NULL)
  }
  
  
  # ----------------------------------------------------------
  # Classes
  # ----------------------------------------------------------
  
  category_lookup <- character(0)
  
  if (
    !is.null(json_data$categories) &&
    length(json_data$categories) > 0
  ) {
    
    ids <- vapply(
      json_data$categories,
      function(x) as.character(x$id),
      character(1)
    )
    
    nomes <- vapply(
      json_data$categories,
      function(x) as.character(x$name),
      character(1)
    )
    
    category_lookup <- setNames(
      nomes,
      ids
    )
  }
  
  
  resultados <- list()
  k <- 0L
  
  
  # ==========================================================
  # Cada imagem registrada no COCO
  # ==========================================================
  
  for (image_info in json_data$images) {
    
    image_coco_id <- as.character(
      image_info$id
    )
    
    largura <- as.integer(
      image_info$width
    )
    
    altura <- as.integer(
      image_info$height
    )
    
    
    # --------------------------------------------------------
    # Descobrir ID espacial
    # --------------------------------------------------------
    
    candidatos <- get_spatial_image_id(
      image_info,
      json_data
    )
    
    encontrados <- candidatos[
      candidatos %in% names(grid_index)
    ]
    
    if (length(encontrados) == 0) {
      
      warning(
        "\nNão encontrei grade para:\n",
        paste(candidatos, collapse = " | "),
        "\nJSON: ", json_file
      )
      
      next
    }
    
    
    spatial_image_id <- encontrados[1]
    
    grid_row <- unname(
      grid_index[spatial_image_id]
    )
    
    grid <- grids[
      grid_row,
    ]
    
    
    # --------------------------------------------------------
    # Selecionar anotações da imagem
    # --------------------------------------------------------
    
    annotations <- Filter(
      function(a) {
        identical(
          as.character(a$image_id),
          image_coco_id
        )
      },
      json_data$annotations
    )
    
    if (length(annotations) == 0) {
      next
    }
    
    
    # ========================================================
    # Cada anotação
    # ========================================================
    
    for (ann in annotations) {
      
      mascara <- decode_coco_rle(
        ann$segmentation
      )
      
      if (
        nrow(mascara) != altura ||
        ncol(mascara) != largura
      ) {
        
        stop(
          "Dimensão da máscara difere de images$width/height em ",
          json_file
        )
      }
      
      
      geom <- mask_to_spatial_polygon(
        mascara = mascara,
        grid_geometry = grid,
        crs = grids,
        crop_to_mask = TRUE
      )
      
      if (is.null(geom)) {
        next
      }
      
      
      # ------------------------------------------------------
      # Opcional:
      # cortar o resultado exatamente pela geometria da grade
      #
      # FALSE é recomendado inicialmente porque suas imagens
      # foram georreferenciadas pelo bbox da grade.
      # ------------------------------------------------------
      
      if (clip_to_grid) {
        
        geom <- suppressWarnings(
          st_intersection(
            geom,
            st_geometry(grid)
          )
        )
        
        if (length(geom) == 0) {
          next
        }
      }
      
      
      # ------------------------------------------------------
      # Metadados
      # ------------------------------------------------------
      
      category_id <- if (!is.null(ann$category_id)) {
        as.character(ann$category_id)
      } else {
        NA_character_
      }
      
      category_name <- if (
        !is.na(category_id) &&
        category_id %in% names(category_lookup)
      ) {
        unname(category_lookup[category_id])
      } else {
        NA_character_
      }
      
      
      score <- if (
        !is.null(ann$score) &&
        length(ann$score) > 0
      ) {
        as.numeric(ann$score)
      } else {
        NA_real_
      }
      
      
      annotation_id <- if (
        !is.null(ann$id)
      ) {
        as.character(ann$id)
      } else {
        NA_character_
      }
      
      
      combio_uid <- if (
        !is.null(ann$combio_uid)
      ) {
        as.character(ann$combio_uid)
      } else {
        NA_character_
      }
      
      
      k <- k + 1L
      
      
      # ------------------------------------------------------
      # Criar feature sf
      # ------------------------------------------------------
      
      resultados[[k]] <- st_sf(
        
        image_id = spatial_image_id,
        
        coco_image_id = image_coco_id,
        
        annotation_id = annotation_id,
        
        category_id = category_id,
        
        category_name = category_name,
        
        score = score,
        
        combio_uid = combio_uid,
        
        area_pixels = sum(mascara),
        
        image_width = largura,
        
        image_height = altura,
        
        # informações da grade
        id_saida = if ("id_saida" %in% names(grid))
          grid$id_saida else NA,
        
        grade_id = if ("grade_id" %in% names(grid))
          grid$grade_id else NA,
        
        NM_MUN = if ("NM_MUN" %in% names(grid))
          as.character(grid$NM_MUN) else NA_character_,
        
        source_json = basename(json_file),
        
        geometry = geom
      )
      
      
      # área espacial
      resultados[[k]]$area_m2 <-
        as.numeric(
          st_area(resultados[[k]])
        )
    }
  }
  
  
  if (length(resultados) == 0) {
    return(NULL)
  }
  
  
  out <- do.call(
    rbind,
    resultados
  )
  
  out
}

