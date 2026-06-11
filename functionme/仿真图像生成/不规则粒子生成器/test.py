class Employee:
    def __init__(self,name,employee_id):
        self.name=name
        self.id=employee_id
    def print_info(self ):
        print(f"该员工姓名为{self.name},工号为{self.id}")
class FullTimeEmployee(Employee):
    def __init__(self,name,employee_id,monthly_salary):
        super().__init__(name,employee_id)
        self.monthly_salary=monthly_salary
    def calculate_monthly_salary(self):
        print(f"全职员工月薪为{self.monthly_salary}")
class PartTimeEmployee(Employee):
    def __init__(self,name,employee_id,daily_salary,work_days):
        super().__init__(name,employee_id)
        self.daily_salary=daily_salary
        self.work_days=work_days
    def calculate_monthly_salary(self):
        print(f"兼职员工月薪为{self.daily_salary*self.work_days}")
worker1=FullTimeEmployee("小张",10086,15000)
worker2=PartTimeEmployee("小李",10087,300,20)
worker1.print_info()
worker1.calculate_monthly_salary()
worker2.print_info()
worker2.calculate_monthly_salary()

