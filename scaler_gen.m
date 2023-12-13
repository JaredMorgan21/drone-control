function numbers = scaler_gen(n)
    persistent counter lastNumbers

    % Initialize counter and lastNumbers if they are empty
    if isempty(counter)
        counter = 0;
        lastNumbers = [];
    end

    % Generate new numbers only every 100 calls
    if counter == 0 || counter >= 10
        while true
            % Generate three random numbers between -1 and 1
            numbers = 2 * rand(3, 1) - 1;

            % Check if their norm is less than or equal to n
            if norm(numbers) <= n
                lastNumbers = numbers; % Store the generated numbers
                counter = 1; % Reset the counter
                break;
            end
        end
    else
        numbers = lastNumbers; % Use the last generated numbers
        counter = counter + 1; % Increment the counter
    end
end