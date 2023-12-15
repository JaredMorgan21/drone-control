function numbers = n_gen(n)
    persistent counter lastNumbers

    % Initialize counter and lastNumbers if they are empty
    if isempty(counter)
        counter = 0;
        lastNumbers = [];
    end

    % Generate new numbers only every 100 calls
    if counter == 0 || counter >= 5000
        while true
            % Generate three random numbers between -1 and 1
            numbers(1) = 0.2 * rand(1) - 0.1;
            numbers(2) = 0.4 * rand(1) - 0.2;
            numbers(3) = 0.02 * rand(1) - 0.01;

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