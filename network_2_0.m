clear all;
% 2.1
p = 0:0.1:0.9;
tN = zeros(length(p), 1);
N = zeros(length(p), 1);
messagesQuantity = 1000;

for i = 1:length(p)
    for j = 1:messagesQuantity
        n = 0;
        receipt = 0;
        while (receipt == 0)
            event = rand();
            if (event > p(i))
                receipt = 1;
            end
            n = n + 1;
        end
        N(i) = N(i) + n;
    end
    N(i) = N(i)/messagesQuantity;
end

for i = 1:length(p)
    tN(i) = 1/(1 - p(i));
end

figure(1);
grid on;
hold on;
plot(p, tN);
plot(p, N);
legend('teorN','N');
hold off;

% 2.2
p = 0:0.1:0.9;
tN1 = zeros(length(p), 1);
N1 = zeros(length(p), 1);
messagesQuantity = 1000;
maxN = 10;

for i = 1:length(p)
    for j = 1:messagesQuantity
        n = 0;
        receipt = 0;
        while (receipt == 0)
            event = rand();
            if ((event > p(i)) || (n == maxN))
                receipt = 1;
            end
            n = n + 1;
        end
        N1(i) = N1(i) + n;
    end
    N1(i) = N1(i)/messagesQuantity;
end

for i = 1:length(p)
    tN1(i) = (1 - p(i)^maxN)/(1 - p(i));
end

figure(2);
grid on;
hold on;
plot(p, tN1);
plot(p, N1);
legend('teorN','N');
hold off;

% 2.3

p = 0:0.1:0.9;
pRe = 0.1;
tN3 = zeros(length(p), length(pRe));
N3 = zeros(length(p), length(pRe));
messagesQuantity = 1000;

for i = 1:length(p)
    for j = 1:length(pRe)
        for k = 1:messagesQuantity
            n = 0;
            receipt = 0;
            while (receipt == 0)
                event = rand();
                if ((event > p(i)))
                    receipt = 1;
                    event = rand();
                    if ((event < pRe(j)))
                        receipt = 0;
                    end
                end
                n = n + 1;
            end
            N3(i,j) = N3(i,j) + n;
        end
        N3(i,j) = N3(i,j)/messagesQuantity;
    end
end

for i = 1:length(p)
    tN3(i,j) = 1/((1 - p(i))*(1 - pRe));
end

figure(3);
grid on;
hold on;
plot(p, tN3);
plot(p, N3);
legend('teorN', 'N');
hold off;

%%%%%%

p = 0:0.1:0.9;
pRe = 0.1;
tN4 = zeros(length(p), length(pRe));
N4 = zeros(length(p), length(pRe));
messagesQuantity = 1000;
maxN = 100;

for i = 1:length(p)
    for j = 1:length(pRe)
        for k = 1:messagesQuantity
            n = 0;
            receipt = 0;
            while (receipt == 0)
                event = rand();
                if ((event > p(i)) || (n == maxN))
                    receipt = 1;
                    event = rand();
                    if ((event < pRe(j)) && (n ~= maxN))
                        receipt = 0;
                    end
                end
                n = n + 1;
            end
            N4(i,j) = N4(i,j) + n;
        end
        N4(i,j) = N4(i,j)/messagesQuantity;
    end
end

for i = 1:length(p)
    tN4(i) = (1 - (1 - (1 - p(i))*(1 - pRe))^maxN)/((1 - p(i))*(1 - pRe));
end

figure(4);
grid on;
hold on;
plot(p, tN4);
plot(p, N4);
legend('teorN', 'N');
hold off;

% 2.4

p = 0.2;
messagesQuantity = 1000;
tau = 2;

t = 1;
sendCounter = 0;
i = 1;
T = 0;

for j = 1:messagesQuantity
    receipt = 0;
    while (receipt == 0)
        spec = 'Ќачало передачи, t = %d, номер сообщени€ = %d\n';
        fprintf(spec, t, j);
        i = i + 1;
        t = t + 1;
        event = rand();
        
        spec = '—ообщение прин€то, t = %d, номер сообщени€ = %d\n';
        fprintf(spec, t, j);
        i = i + 1;
        t = t + tau - 1;
        if (event > p)
            receipt = 1;
            spec = 'ѕолучена положительна€ квитанци€, t = %d, номер сообщени€ = %d\n\n';
            fprintf(spec, t, j);
            i = i + 1;
        else    
            spec = 'ѕолучена отрицательна€ квитанци€, t = %d, номер сообщени€ = %d\n\n';
            fprintf(spec, t, j);
            i = i + 1;
        end
        t = t + 1;
        sendCounter = sendCounter + 1;
    end
end

tN = (1 - p)/(1 + tau);
disp(tN);

N = messagesQuantity/t;
disp(N);


% 2.5

p = 0.2;
messagesQuantity = 1000;
tau = 2;

t = 2;
sendCounter = 0;
sendedMessages = 0;
i = 1;
messages = ones(messagesQuantity, 1);
receipts = ones(messagesQuantity, 1);
sendMemory = zeros(1, tau);
receiptMemory = zeros(1, tau);

while (sendedMessages < messagesQuantity)
    
    % отправить сообщение
    for i = 1:length(messages)
        if (messages(i) == 1)
            spec = 'Ќачало передачи, t = %d, номер сообщени€ = %d\n';
            fprintf(spec, t, i);
            event = rand();
            if (event > p)
                messages(i) = t + 1;
            else
                messages(i) = -(t + 1);
            end
            break;
        end
    end
    
    % прин€ть сообщение
    for i = 1:length(messages)
        flag = 1;
        for k = 1:length(sendMemory)
            if (i == sendMemory(k))
                spec = 'ѕропуск сообщени€ после ошибки, t = %d, номер сообщени€ = %d\n';
                fprintf(spec, t, i);
                flag = 0;
                sendMemory(k) = 0;
                break;
            end
        end
        if (flag == 1)
            if ((abs(messages(i)) == t) && (abs(messages(i)) ~= 1))
                if (messages(i) > 0)
                    receipts(i) = messages(i) - 1 + tau;
                    spec = 'Cообщение прин€то, t = %d, номер сообщени€ = %d\n';
                    fprintf(spec, t, i);
                else
                    receipts(i) = messages(i) + 1 - tau;
                    spec = 'Cообщение прин€то с ошибкой, t = %d, номер сообщени€ = %d\n';
                    fprintf(spec, t, i);
                    for b = 1:tau
                        sendMemory(b) = i + b;
                        receiptMemory(b) = i + b;
                    end
                end
                break;
            end
        else
            break;
        end
    end
    
    % прин€ть квитанцию
    for i = 1:length(receipts)
        flag = 1;
            for k = 1:length(receiptMemory)
                if (i == receiptMemory(k))
                    spec = 'ѕолучена отрицательна€ квитанци€ пропуска сообщени€, t = %d, номер сообщени€ = %d\n';
                    fprintf(spec, t, i);
                    flag = 0;
                    receiptMemory(k) = 0;
                    break;
                end
            end
        if (flag == 1)
            if ((abs(receipts(i)) == t) && (abs(receipts(i)) ~= 1))
                if (receipts(i) > 0)
                    receipts(i) = 0;
                    sendedMessages = sendedMessages + 1;
                    messages(i) = 0;
                    spec = 'ѕолучена положительна€ квитанци€, t = %d, номер сообщени€ = %d\n';
                    fprintf(spec, t, i);
                else
                    spec = 'ѕолучена отрицательна€ квитанци€, t = %d, номер сообщени€ = %d\n';
                    fprintf(spec, t, i);
                    for j = i:length(messages)
                        messages(j) = 1;
                    end
                    for j = i:length(receipts)
                        receipts(j) = 1;
                    end
                end
                break;
            end
        else
            break;
        end
    end
    fprintf('\n');
    t = t + 1;
    
end

tN = (1 - p)/(1 + tau*p);
disp(tN);

N = messagesQuantity/t;
disp(N);








