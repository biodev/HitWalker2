#' Neo4j Interface Functions
#'
#' These functions utilities which assist with working with a Neo4j database
#'
#' @name neo4j_helpers
NULL

make.property.str.node <- function(prop.names, dta, pos, array.delim)
    {
        base.str <- paste0("{name:line[",pos,"]")
        
        for(i in prop.names)
        {
            use.name <- sapply(strsplit(i, "\\."), function(x) x[-1])
            use.pos <- which(names(dta) == i)-1
            
            stopifnot(length(use.pos) == 1)
            
            if (is.integer(dta[,i]))
            {
                base.str <- paste0(base.str, ",",use.name,":toInt(line[",use.pos,"])")
            }else if (is.numeric(dta[,i]))
            {
                base.str <- paste0(base.str, ",",use.name,":toFloat(line[",use.pos,"])")
            }else{
                
                if (is.null(array.delim) == F && is.character(array.delim) == T && length(array.delim) == 1)
                {
                    if (any(grepl(array.delim, dta[,i])))
                    {
                        base.str <- paste0(base.str, ",",use.name,":split(line[",use.pos,"], '",array.delim,"')")
                    }else{
                        base.str <- paste0(base.str, ",",use.name,":line[",use.pos,"]")
                    }
                }else{
                    base.str <- paste0(base.str, ",",use.name,":line[",use.pos,"]")
                }
            }
        }
        
        base.str <- paste0(base.str, "}")
        
        return(base.str)
    }

get.or.create.constraints <- function(neo.path=NULL, node.names)
{
    
    if (missing(neo.path) == F && is.null(neo.path) == F && all(is.na(neo.path)) == F)
    {
       use.neo.path <- file.path(neo.path, "bin", "neo4j-shell")
    }else{
        
        use.neo.path <- "neo4j-shell"
    }
    
    #before executing, things will go A LOT faster if constraints/indexes are available, check first and if not, then create them
    
    found.schema <- system(paste(use.neo.path, "-c schema"), intern=T)
    
    constraint.section <- grep("Constraints", found.schema)
    
    if (length(constraint.section) > 0)
    {
        found.constraints <- found.schema[(constraint.section+1):length(found.schema)]
    
        present.constraints <- sapply(node.names, function(x)
               {
                    any(grepl(paste0("\\s+ON\\s+\\(",tolower(x),":",x,"\\) ASSERT\\s+", tolower(x), ".name IS UNIQUE"), found.constraints))
               })
    }else{
        #no constraints were found so set both notes to present.constraints=F
        present.constraints <- rep(F,2)
        names(present.constraints) <- node.names
    }
    
    if (all(present.constraints))
    {
        message("Found all necessary CONSTRAINTS, continuing...")
        missing.constraints <- character(0)
    }else{
        missing.constraints <- names(present.constraints)[present.constraints == F]
        message(paste("Missing CONSTRAINT(s) for:", paste(missing.constraints, collapse=","), "will build them first and continue with loading..."))
    }
    
    if (length(missing.constraints) > 0)
    {
        constr.str <- paste0(use.neo.path, " -c 'CREATE CONSTRAINT ON (",tolower(missing.constraints),":",missing.constraints,") ASSERT ",tolower(missing.constraints),".name IS UNIQUE;'")
    }else{
        constr.str <- character(0)
    }
    
    return (constr.str)
}

#' @rdname neo4j_helpers
#' @param .data A \code{data.frame} containing data to load where the first two columns are nodes and each row is an implied edge. The first two
#' column names are taken to be the names of the nodes.  The following columns should either be of the form 'node.val' or simply 'val'.  The former
#' provides a property for the specified node, while the latter specifies an edge property.
#' @param edge.name Name of the implied edge
#' @param commit.size Number of relationships to commit in a single transaction when loading.
#' @param dry.run If TRUE, then the actual statments will be printed to the screen but not executed.
#' @param array.delim If a property should encode a Neo4j array, this specifies how the values should be delimited.
#' @param unique.rels If TRUE, relationships will be required to be unique.
#' @param merge.from Specifies whether the first node should be merged, that is added if it doesn't exist, or simply matched.
#' @param merge.to Species whether the second node should be merged similar to merge.from.
#' @return Nothing, as a side effect the specified Neo4j database is populated.
load.neo4j <- function(.data, edge.name=NULL, commit.size=1000, neo.path=NULL, dry.run=F, array.delim="&", unique.rels=T, merge.from=T, merge.to=T)
{
    if (missing(edge.name) || is.null(edge.name))
    {
        load.neo4j.node(.data=.data, commit.size=commit.size, neo.path=neo.path, dry.run=dry.run, array.delim=array.delim)
        
    }else{
        load.neo4j.edge(.data=.data, edge.name=edge.name, commit.size=commit.size, neo.path=neo.path, dry.run=dry.run, array.delim=array.delim, unique.rels=unique.rels, merge.from=merge.from, merge.to=merge.to)
    }
}

load.neo4j.node <- function(.data, commit.size=1000, neo.path=NULL, dry.run=F, array.delim="&")
{
    if (any(is.na(.data)))
    {
        stop("ERROR: .data currently cannot have NAs, please remove and try again...")
    }
    
    base.node.names <- names(.data)[1]
    
    node.names <- capwords(base.node.names)
    
    node.props.dta <- lapply(base.node.names, function(x)
                         {
                            return(names(.data)[grep(paste0(x,"\\."), names(.data))])
                         })
    
    use.temp <- tempfile()
    
    write.table(.data, sep="\t", file=use.temp, col.names=F, row.names=F, quote=F)
    
    cypher.temp <- tempfile()
    
    constr.str <- get.or.create.constraints(neo.path, node.names)
    
    cypher.stats <- c(
        paste0("USING PERIODIC COMMIT ",commit.size),
        paste0("LOAD CSV FROM 'file://", use.temp,"' AS line FIELDTERMINATOR '\t'"),
        paste0("MERGE (n:",node.names[1], make.property.str.node(node.props.dta[[1]], .data, 0, array.delim),");")
    )
    
    writeLines(cypher.stats, con=cypher.temp)
    
    if (missing(neo.path) == F && is.null(neo.path) == F && all(is.na(neo.path)) == F)
    {
       use.neo.path <- file.path(neo.path, "bin", "neo4j-shell")
    }else{
        
        use.neo.path <- "neo4j-shell"
    }
    
    if (dry.run == F)
    {
        for (i in constr.str)
        {
            system(i)
        }
        
        message("Loading into Neo4j...")
        system(paste(use.neo.path, "-file", cypher.temp))
    }else{
        
        for (i in constr.str)
        {
            message(i)
        }
        
        message(paste(use.neo.path, "-file", cypher.temp))
    }
}

#data.frame should be in the form:
#data.frame(node1, node2, node[12].property, edge property(no x'.'y just name))
load.neo4j.edge <- function(.data, edge.name=NULL, commit.size=1000, neo.path=NULL, dry.run=F, array.delim="&", unique.rels=T, merge.from=T, merge.to=T)
{
    message("Preprocessing...")
    
    if (missing(edge.name) || is.null(edge.name) || is.na(edge.name))
    {
        stop("ERROR: Need to supply an edge name")
    }else if (is.character(edge.name) == F || length(edge.name) != 1)
    {
        stop("ERROR: edge.name needs to be a single string value")
    }
    
    if (ncol(.data) < 2 || is.null(names(.data)))
    {
        stop("ERROR: .data needs to have at least two columns and be named")
    }
    
    if (any(is.na(.data)))
    {
        stop("ERROR: .data currently cannot have NAs, please remove and try again...")
    }
    
    
    .make.property.str.edge <- function(prop.names, dta, array.delim)
    {
        if (length(prop.names) == 0)
        {
            return("")
        }
        
        base.str <- "{"
        
        for(i in prop.names)
        {
            if (which(prop.names == i) > 1)
            {
                base.str <- paste0(base.str, ",")
            }
            
            use.pos <- which(names(dta) == i)-1
            
            stopifnot(length(use.pos) == 1)
            
            if (is.integer(dta[,i]))
            {
                base.str <- paste0(base.str,i,":toInt(line[",use.pos,"])")
            }else if (is.numeric(dta[,i]))
            {
                base.str <- paste0(base.str,i,":toFloat(line[",use.pos,"])")
            }else{
                
                if (is.null(array.delim) == F && is.character(array.delim) == T && length(array.delim) == 1)
                {
                    if (any(grepl(array.delim, dta[,i])))
                    {   
                        base.str <- paste0(base.str,i,":split(line[",use.pos,"], '", array.delim,"')")
                    }else{
                        base.str <- paste0(base.str,i,":line[",use.pos,"]")
                    }
                }else{
                
                    base.str <- paste0(base.str,i,":line[",use.pos,"]")
                
                }
            }
        }
        
        base.str <- paste0(base.str, "}")
        
        return(base.str)
    }
    
    base.node.names <- names(.data)[1:2]
    
    node.names <- capwords(base.node.names)
    
    node.props.dta <- lapply(base.node.names, function(x)
                         {
                            return(names(.data)[grep(paste0(x,"\\."), names(.data))])
                         })
    
    edge.props <- setdiff(names(.data), c(base.node.names, unlist(node.props.dta)))
    
    use.temp <- tempfile()
    
    write.table(.data, sep="\t", file=use.temp, col.names=F, row.names=F, quote=F)
    
    cypher.temp <- tempfile()
    
    constr.str <- get.or.create.constraints(neo.path, node.names)
    
    rel.unique.str <- ifelse(unique.rels, "UNIQUE", "")
    
    if (merge.from){
        from_str = "MERGE"
    }else{
        from_str = "MATCH"
    }
    
    if (merge.to){
        to_str = "MERGE"
    }else{
        to_str = "MATCH"
    }
    
    cypher.stats <- c(
        paste0("USING PERIODIC COMMIT ",commit.size),
        paste0("LOAD CSV FROM 'file://", use.temp,"' AS line FIELDTERMINATOR '\t'"),
        paste0(from_str, " (n:",node.names[1], make.property.str.node(node.props.dta[[1]], .data, 0, array.delim),")"),
        paste0(to_str, " (m:",node.names[2],make.property.str.node(node.props.dta[[2]], .data, 1, array.delim),")"),
        paste0("CREATE ",rel.unique.str," (n)-[:",edge.name,.make.property.str.edge(edge.props, .data, array.delim),"]->(m);")
    )
    
    writeLines(cypher.stats, con=cypher.temp)
    
    if (missing(neo.path) == F && is.null(neo.path) == F && all(is.na(neo.path)) == F)
    {
       use.neo.path <- file.path(neo.path, "bin", "neo4j-shell")
    }else{
        
        use.neo.path <- "neo4j-shell"
    }
    
    if (dry.run == F)
    {
        for (i in constr.str)
        {
            system(i)
        }
        
        message("Loading into Neo4j...")
        system(paste(use.neo.path, "-file", cypher.temp))
    }else{
        
        for (i in constr.str)
        {
            message(i)
        }
        
        message(paste(use.neo.path, "-file", cypher.temp))
    }
}

clean.neo4j.res <- function(result)
{
    if (length(result) == 6)
    {
        return(NULL)
    }else{
        
        use.res <- result[seq(from=4, to=length(result)-3)]
    
        #also get the header
        use.res <- append(result[2], use.res)
        
        use.res <- gsub("\"", "", gsub("|\\[|\\]|\\s", "", use.res), fixed=T)
        
        split.res <- sapply(strsplit(use.res, "\\|"), function(x) x[x!=""])
        
        if (class(split.res) == "character")
        {
            return(split.res[-1])
        }else if (class(split.res) == "matrix")
        {
            header <- split.res[,1]
            split.res <- split.res[,-1]
            ret.dta <- data.frame(t(split.res), stringsAsFactors=F)
            names(ret.dta) <- make.names(header)
            return(ret.dta)
        }else{
            stop("ERROR: Unexpected type for result")   
        }
    }
    
    
}

#' @rdname neo4j_helpers
#' @param neo.path If \code{neo.path} is specified, the neo4j-shell executable is expected at neo.path/bin/neo4j-shell.  Otherwise it is expected to be part of your path.
#' @return A \code{igraph} object depicting the database.
compute.graph.structure <- function(neo.path=NULL)
{
    message("Getting labels from DB")
    
    if (missing(neo.path) == F && is.null(neo.path) == F && all(is.na(neo.path)) == F)
    {
       use.neo.path <- file.path(neo.path, "bin", "neo4j-shell")
    }else{
        
        use.neo.path <- "neo4j-shell"
    }
    
    label.query <- "'MATCH (n) RETURN DISTINCT LABELS(n);'"
    
    labs <- system(paste(use.neo.path, "-c", label.query), intern=T)
    
    use.labs <- clean.neo4j.res(labs)
    
    lab.list <- lapply(use.labs, function(i)
           {
                message(paste("Starting", i))
                cur.query <- paste0('MATCH (n:', i ,')-[r]->(m) RETURN DISTINCT LABELS(n) AS from_node, LABELS(m) AS to_node, TYPE(r) AS type;')
        
                cur.res <- system(paste0(use.neo.path, " -c ","'",cur.query, "'"), intern=T)
                
                return(clean.neo4j.res(cur.res))
           })
    
    lab.dta <- do.call("rbind", lab.list)
    
    #make an igraph object
    
    return(graph.data.frame(lab.dta))
}

#neo.path=normalizePath("neo4j-community-2.1.6/")
make.graph.struct <- function(neo.graph, graph.struct.path="test_graph_struct.json")
{
    
    if (missing(graph.struct.path) || is.null(graph.struct.path) || all(is.na(graph.struct.path)) || (is.character(graph.struct.path) == F))
    {
       stop("ERROR: Need to supply a valid path for the graph_struct file")
    }
    
    if (missing(neo.graph) || is.null(neo.graph) || all(is.na(neo.graph)) || (class(neo.graph) != "igraph"))
    {
        stop("ERROR: neo.graph needs to be an igraph object")
    }
    
    lab.graph <- neo.graph
    
    un.neo.graph <- as.undirected(neo.graph, edge.attr.comb="first")

    graph.list <- lapply(get.adjedgelist(un.neo.graph), function(x)
                         {
                            unique.x <- unique(x)
                            
                            temp.list <- lapply(unique.x, function(y)
                                   {
                                        return(V(un.neo.graph)[inc(E(un.neo.graph)[y])]$name)
                                   })
                            
                            names(temp.list) <- E(un.neo.graph)[unique.x]$type
                            return(temp.list)
                         })
    
    pretty.graph.list <- mapply(function(x,y){
        
        lapply(x, function(z)
               {
                    if (length(z) > 1)
                    {
                        return(z[z != y])
                    }else{
                        return(z)
                    }
                    
               })
    }, graph.list, names(graph.list))
    
    write(toJSON(pretty.graph.list), file=graph.struct.path)
    
}

read.graph_struct  <- function(graph_struct="/var/www/hitwalker_2_inst/graph_struct.json")
{
    
    json.obj <- fromJSON(file=graph_struct)
    
    #warnings are about the rownames...
    json.dta <- suppressWarnings(data.frame(do.call(rbind, lapply(names(json.obj), function(x) {
        if (length(json.obj[[x]]) > 0)
        {
            temp.dta <- cbind(x, names(json.obj[[x]]), unlist(json.obj[[x]]))
        }else{
            return(NULL)
        }
        
    })), stringsAsFactors=F))
    names(json.dta) <- c("from", "label", "to")
    json.dta <- json.dta[!duplicated(json.dta$label),]
    json.graph <- graph.data.frame(json.dta[,c("from", "to", "label")], directed=F)
    return(json.graph)
    #plot(json.graph)
}