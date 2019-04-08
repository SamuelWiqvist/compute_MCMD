function game(grid::Array, wrap_around::Bool)

    current_place = 1

    run_game = true
    nbr_tosses = 0

    grid[current_place] = grid[current_place]+1

    while run_game

        nbr_tosses = nbr_tosses + 1
        steps = rand(1:6)
        current_place = current_place+steps

        if current_place == 81
            run_game = false
            grid[current_place] = grid[current_place]+1
        elseif current_place > 81
            if wrap_around # wrap around
                current_place = current_place-81
            else
                current_place = 81-(current_place-81)
            end
            grid[current_place] = grid[current_place] + 1
        else
            grid[current_place] = grid[current_place]+1
        end

    end

    return nbr_tosses
end
