
macro facet_df(df, args...)
    if length(args) < 2
        error("@facet_df requires at least one grouping column and a plot expression.")
    end

    # Extract column symbols (handles :col as QuoteNode and col as Symbol)
    cols = Symbol[arg isa QuoteNode ? arg.value : arg for arg in args[1:(end-1)]]
    expr = args[end]

    # Escape the entire generated block at the top level
    return esc(quote
        [
            @df subdf $expr
            for subdf in DataFrames.groupby(DataFrames.sort($df, $cols), $cols)
        ]
    end)
end
